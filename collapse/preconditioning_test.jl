#################################################################
#                   Includes and Modules                        #
#################################################################

using LinearAlgebra, Arpack, Random, SparseArrays
using BenchmarkTools, Printf, TimerOutputs
Random.seed!(12345)

#include("../m1Accel.jl")

include("../HOLaGraf_lsmr.jl")
using .HOLaGraf_lsmr
include("../generateDelauney.jl")

# Plotting thingies 
using Plots, ColorSchemes, Plots.PlotMeasures,  LaTeXStrings
pgfplotsx()
theme(:mute)
Plots.scalefontsizes(1.25)
cols=ColorSchemes.Spectral_11;

# Sampling
using StatsBase
using GR: delaunay

#################################################################
#                End of Includes and Modules                    #
#################################################################


#################################################################
#                   Functions of Sets                          #
#################################################################

"""
    getEdges2Trig( edges2, trigs2 )

Create a vector of sets of triangles adjacent to edges
"""
function getEdges2Trig( edges2, trigs2 )
      edg2Trig = [];
      for i in axes(edges2, 1)
            tmp = [];
            for j in axes(trigs2, 1)
                  if ((trigs2[j, 1] == edges2[i, 1]) && (trigs2[j, 2] == edges2[i, 2])) || ((trigs2[j, 1] == edges2[i, 1]) && (trigs2[j, 3] == edges2[i, 2])) || ((trigs2[j, 2] == edges2[i, 1]) && (trigs2[j, 3] == edges2[i, 2]))
                        tmp=[tmp; j];
                  end
            end
            edg2Trig = [edg2Trig; Set(tmp)];
      end
      return edg2Trig
end
      
"""
    getTrig2Edg( edges2, trigs2, edg2Trig )

create a vector of sets of edges adjacent to triangles
"""
function getTrig2Edg( edges2, trigs2, edg2Trig )
      trig2Edg = [ Set{Int64}([]) for i in axes(trigs2, 1) ];
      for i in axes(edges2, 1)
            for j in edg2Trig[i]
                  push!(trig2Edg[j], i)
            end
      end
      return trig2Edg
end
      
#################################################################
#             END of Functions of Sets                          #
#################################################################



#################################################################
#                Functions of Delaunay                          #
#################################################################

function sparseDelaunay( ; N = 10, ν = 0.4 )
      n  = N + 4;

      points, edges, trigs = generateDelauney( N )
      ν_init = size(edges, 1 ) / ( ( N + 4) * ( N + 3 ) / 2 );

      backlash = Int( - size( edges, 1 ) + round( ν * n * (n-1) / 2 ) ) 

      if backlash < 0
            for i in 1 : -backlash
                  ind = getIndx2Kill( edges ) ;
                  edges, trigs = killEdge(ind, n, edges, trigs)
            end
      else
            allEdges = Array{Integer}(undef, 0, 2)
            for i in 1:(n-1) 
                  for j in (i+1):n
                        allEdges = [ allEdges; i j  ]
                  end
            end
            for i in 1:size(edges, 1)
                  indx = findall(all(allEdges .== edges[i, :]', dims=2))[1][1]
                  allEdges = allEdges[ 1:size(allEdges, 1) .!= indx, : ]
            end
            for i in 1 : backlash      
                  ind, allEdges = getNewEdge2(n, edges, allEdges);
                  edges, trigs = addEdge(ind, n, edges, trigs)
            end
      end

      return n, points, ν_init, edges, trigs

end

#################################################################
#            END of Functions of Delaunay                       #
#################################################################


function greedyCollapse( edges, trigs )
      edge2Trig = getEdges2Trig(edges, trigs)
      trig2Edg = getTrig2Edg(edges, trigs, edge2Trig)

      Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
      Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

      Σ = [];
      Τ = [ ]; 
      while !isempty(Free)
            σ = pop!(Free) 
            Σ = [ Σ; σ ]
            τ = first( edge2Trig[ σ ] ) 
            Τ = [ Τ ; τ ]
            for ξ in trig2Edg[ τ ]
                  setdiff!( edge2Trig[ ξ ], Set([ τ ]) )
                  Ls[ ξ ] = Ls[ ξ ] - 1
                  ( Ls[ ξ ] == 1 ) && ( union!(Free, Set( [ ξ ] ) ) )
                  ( Ls[ ξ ] == 0 ) && ( setdiff!(Free, Set( [ ξ ] ) ) )
            end
      end
      return sum( Ls ) == 0, Σ, Τ, edge2Trig, trig2Edg, Ls, Free  
end


function generateComplete( n )
      edges = Array{Integer}(undef, 0, 2)
      for i in 1 : n-1
            for j in i+1 : n
                  edges = [edges; i j ]
            end
      end
      trigs = Array{Integer}(undef, 0, 3)
      for i in 1 : n-2
            for j in i+1 : n-1
                  for k in j+1 : n
                        trigs = [ trigs; i j k ]
                  end
            end
      end

      return n, edges, trigs
end

function indicator( ind, m )
      res = zeros(m)
      res[ind] .= 1
      return diagm(res)
end



function B2fromTrig(edges, trigs)
    
      m = size(edges, 1);
      if length(trigs)==1
          return spzeros(m, 1)
      end
      del = size(trigs, 1);
      B2 = spzeros(m, del);
      
      for i in 1:del
          B2[findfirst( all( [trigs[i, 1], trigs[i, 2]]' .== edges, dims=2 )[:, 1] ), i]=1;
          B2[findfirst( all( [trigs[i, 1], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=-1;
          B2[findfirst( all( [trigs[i, 2], trigs[i, 3]]' .== edges, dims=2 )[:, 1] ), i]=1;
      end 
      return B2
end
 
function Perm(v)
      m = size(v, 1)
      res = zeros( m, m )
      for i in 1 : m 
            res[ v[i], i ] = 1
      end
      return res
end

function condPlus( A; thr = 1e-8 )
      m = size(A, 1);
      #σ = svds(L1up, nsv = m - 1)[1].S;
      σ = svd(Matrix(A)).S;
      return maximum( σ ) / minimum( σ[ abs.(σ) .> thr ])
end



n, edges, trigs = generateComplete(4)
w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25
W = diagm(vec(sqrt.(w)))
m = size( edges, 1 )
Δ = size( trigs, 1 )

B2 = B2fromTrig( edges, trigs )
Lu =  B2 * W * W * B2'

original = condPlus( Lu )

filter = [1; 3; 4]
Π = indicator( filter, Δ )
trigs2 = trigs[filter, :]
edge2Trig2 = getEdges2Trig( edges, trigs2 )
trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )

fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2
D = B2 * W * Π
E = B2 * W




condPlus( Lu ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') ), 
condPlus( P2' * pinv(D) *  Lu * pinv(D') * P2 ), condPlus( pinv(D) *  Lu * pinv(D') ),
condPlus( Π * pinv(E) *  Lu * pinv(E') * Π ), condPlus( Π * pinv(E) * E * E' * pinv(E') * Π )

condPlus( Π * pinv(E) * E )
condPlus( E )

condPlus( pinv(E) * E )



X = rand(9, 15)
Π = indicator( sample(1:15, 10), 15 )
condPlus(X), condPlus( Π * pinv(X) * X )

Y = [X; X[1,:]']
condPlus(Y), condPlus( Π * pinv(Y) * Y )


function singulars( X )
      return eigvals( Matrix(X' * X) )
end


U, S, V = svd(E)
S = diagm(S)



singulars( E )
[ singulars( S ) 1 ./singulars( pinv(S) ) ]
singulars( Π * pinv( E ) * E )

condPlus( Π * pinv(E) * E )
condPlus( Π * V * pinv(S) * S )
condPlus( E )

condPlus( pinv(E) )
condPlus( Π * V )

condPlus( S )
condPlus( pinv(S) * S )




n = 6
edges = [ 1 2; 1 3; 1 4; 1 5; 2 3; 2 4; 2 6; 3 5; 3 6; 4 5; 4 6; 5 6]
trigs = [ 1 2 3; 1 2 4; 1 3 5; 1 4 5; 2 3 6; 2 4 6; 3 5 6; 4 5 6]

w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25
W = diagm(vec(sqrt.(w)))
m = size( edges, 1 )
Δ = size( trigs, 1 )

B2 = B2fromTrig( edges, trigs )
Lu =  B2 * W * W * B2'

condPlus( Lu )

ind = sample( 1 : size(trigs, 1) )
filter = findall( 1 : size(trigs, 1) .!= ind)

Π = indicator( filter, Δ )
trigs2 = trigs[filter, :]
edge2Trig2 = getEdges2Trig( edges, trigs2 )
trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )

fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2
D = B2 * W * Π
E = B2 * W

fl, condPlus( Lu ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') )

U, S, V = svd(E)
S = diagm(S)

X = V' * Π * V
SI = pinv(S) * S










n, edges, trigs = generateComplete(10)
w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25
W = diagm(vec(sqrt.(w)))
m = size( edges, 1 )
Δ = size( trigs, 1 )

B2 = B2fromTrig( edges, trigs )

edge2Trig = getEdges2Trig( edges, trigs )
trig2Edge = getTrig2Edg( edges, trigs, edge2Trig )

num = Int( round( 2/3 * m ) ) - 3

while true
      filter = sort( sample( 1:m, num ) )
      Π = indicator( filter, Δ )

      trigs2 = trigs[filter, :]
      edge2Trig2 = getEdges2Trig( edges, trigs2 )
      trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )
      fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )
      if fl
            break
      end
end

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2[:, 1:m]
Lu = B2 * W * W * B2'


condPlus( Lu ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') )
norm( Lu - I(m) )/(m^2), norm( pinv(C) * P1 * Lu * P1' * pinv(C') - I(m) )/(m^2)

norm( Lu - I(m) ), norm( pinv(C) * P1 * Lu * P1' * pinv(C') - I(Δ) )




















n, edges, trigs = generateComplete(4)
w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25
W = diagm(vec(sqrt.(w)))
m = size( edges, 1 )
Δ = size( trigs, 1 )

B2 = B2fromTrig( edges, trigs )
Lu =  B2 * W * W * B2'

original = condPlus( Lu )

filter = [1; 2; 3]
Π = indicator( filter, Δ )
trigs2 = trigs[filter, :]
edge2Trig2 = getEdges2Trig( edges, trigs2 )
trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )

fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2
D = B2 * W * Π
E = B2 * W




w, condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') )/condPlus( Lu ) 























n = 7
edges = [ 1 2; 1 3; 1 4; 1 6; 2 3; 2 5; 2 6; 3 4; 3 5; 4 5; 4 6; 4 7; 5 6; 5 7; 6 7]
trigs = [ 1 2 3; 1 2 6; 1 3 4; 1 4 6; 2 3 5; 2 5 6; 3 4 5; 4 5 6; 4 5 7; 4 6 7; 5 6 7]
w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25
W = diagm(vec(sqrt.(w)))
m = size( edges, 1 )
Δ = size( trigs, 1 )

B2 = B2fromTrig( edges, trigs )
Lu =  B2 * W * W * B2'

original = condPlus( Lu )


ind = sample(1:8)
filter = findall( 1:Δ .!= ind )
ind2 = sample(9:11)
filter = filter[findall( filter .!= ind2 )]

#filter = [ 1; 2; 3; 5; 6; 7; 8; 10; 11 ]
filter = [ 2; 3; 5; 6; 7; 10; 11]
Π = indicator( filter, Δ )

trigs2 = trigs[filter, :]
edge2Trig2 = getEdges2Trig( edges, trigs2 )
trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )

fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )

Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

P1 = Perm( Σ_full )
P2 = Perm( Τ_full )

C = P1 * B2 * W * Π * P2
D = B2 * W * Π
E = B2 * W

fl, w, condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') )/condPlus( Lu ) 
w[ ind ] / minimum( w[ 1 : 8 ] ), w[ ind2 ] / minimum( w[ 9 : 11 ] ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') )/condPlus( Lu ) 































