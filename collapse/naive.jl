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

function greedyCollapseShort( edge2Trig, trig2Edg, perm, edges)
      
      Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
      Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

      Σ = [];
      Τ = [ ]; 
      while !isempty(Free)
            σ = pop!(Free) 
            Σ = [ Σ; σ ]
            τ = first( edge2Trig[ σ ] ) 
            Τ = [ Τ ; τ ]
            in = findall(perm .== τ)[1][1] 
            for ξ in trig2Edg[ in ]
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

function condPlus( A; thr = 1e-4 )
      m = size(A, 1);
      #σ = svds(L1up, nsv = m - 1)[1].S;
      σ = svd(Matrix(A)).S;
      return maximum( σ ) / minimum( σ[ abs.(σ) .> thr ])
end





n, edges, trigs = generateComplete(8)

#n = 7
#edges = [ 1 2; 1 3; 1 4; 1 6; 2 3; 2 5; 2 6; 3 4; 3 5; 4 5; 4 6; 4 7; 5 6; 5 7; 6 7]
#trigs = [ 1 2 3; 1 2 6; 1 3 4; 1 4 6; 2 3 5; 2 5 6; 3 4 5; 4 5 6; 4 5 7; 4 6 7; 5 6 7]

#w = 0.99*rand( size(trigs, 1), 1 ) .+ 0.0001
#w = abs.( randn( size(trigs, 1)) ) / 3 .+ 0.0001
w = abs.( rand( Cauchy(0, 1) , size(trigs, 1) ) )
#w[ w.> 3 ] .= 3

#w =[ 1.0; 0.01; 0.01; 1.0]
#w = [ 0.5075182401719001; 0.875870266151905; 0.30441538509953425; 0.9924611971080138; 0.7477281699815993; 0.4376070266421683; 0.46207272132020705; 0.6159243377415332; 0.5775446077858322; 0.6641950963632279; 0.7617574896668902]

#perm = sortperm( w, rev = true )

m = size( edges, 1 )
Δ = size( trigs, 1 )

#perm = sortperm( w, dims = 1, rev = true )
perm = sortperm( w, rev = true )
W = diagm(vec(sqrt.(w)))

B2 = B2fromTrig( edges, trigs )
Lu =  B2 * W * W * B2'

original = condPlus( Lu )

edg2Trig = getEdges2Trig( edges, trigs )
trig2Edge = getTrig2Edg( edges, trigs, edg2Trig )




sub = [ ] 
ind = 1

subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
subTrig2Edg = [ ]

while ind <= Δ
      global subEdg2Trigs, subTrig2Edg, trig2Edg, sub, ind
      tmpSubEdg2Trigs = deepcopy( subEdg2Trigs )
      tmpSubTrig2Edg = deepcopy( subTrig2Edg )
      tmpSub = deepcopy( sub )

      tmpSub = [ tmpSub; perm[ind] ]
      tmpSubTrig2Edg = [ tmpSubTrig2Edg; trig2Edge[ perm[ind] ] ]
      for e in trig2Edge[ perm[ind] ]
            tmpSubEdg2Trigs[e] = union( tmpSubEdg2Trigs[e], perm[ind] )
      end

      fl, Σ, Τ, edge2Trig, trig2Edg, Ls, Free = greedyCollapseShort( tmpSubEdg2Trigs, tmpSubTrig2Edg, tmpSub )
      
      if fl
            sub = [ sub; perm[ind] ]
            subTrig2Edg = [ subTrig2Edg; trig2Edge[ perm[ind] ] ]
            for e in trig2Edge[ perm[ind] ]
                  subEdg2Trigs[e] = union( subEdg2Trigs[e], perm[ind] )
            end
      end

      ind = ind + 1
end




filter = sub
Π = indicator( sort(filter), Δ )
trigs2 = trigs[ sort(filter), :]
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

fl, condPlus( Lu ), condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') ), norm(sqrt.(W)*Π)/norm(sqrt.(W)) #size(nullspace( Π ), 2), Δ- size(nullspace( D ), 2), minimum(w) 


sqrt.( abs.( singulars( Lu ) ) )
sqrt.( abs.( singulars( pinv(C) * P1 * Lu * P1' * pinv(C') ) ) )


bitmask(z) = digits(z, base=2, pad=20) #|> reverse;


full = 2^Δ
flags = zeros( full, 1 )
condNums = zeros( full, 1 )
kers = zeros( full, 1 )

dimTest = size(nullspace(W*B2'), 2)

for i in 2:full

      global full, flags, condNums, Δ, edges, Lu, W, B2
      
      filter = findall( bitmask(i) .== 1 )
      Π = indicator( sort(filter), Δ )
      trigs2 = trigs[ sort(filter), :]
      edge2Trig2 = getEdges2Trig( edges, trigs2 )
      trig2Edge2 = getTrig2Edg( edges, trigs2, edge2Trig2 )

      fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges, trigs2 )

      Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
      Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

      P1 = Perm( Σ_full )
      P2 = Perm( Τ_full )

      C = P1 * B2 * W * Π * P2
      
      flags[i] = fl
      condNums[i] = condPlus( pinv(C) * P1 * Lu * P1' * pinv(C') ) 
      kers[i] = ( size(nullspace(Π*W*B2'), 2) == size(nullspace(W*B2'), 2) )
end

















tol = 1e-3

A = Lu
b = Lu * ones(28)



x0 = zeros( size(A, 2), 1 )
      
r0 = b - A * x0
p0 = A' * r0
s0 = p0

γ = norm(s0)^2
β = 2.0


qi = A * p0
ai = γ / norm(qi)^2
x0 = x0 + ai * p0
r0 = r0 - ai * qi

abs( β - 1 ) < tol

s0 = A' * r0
β = norm( s0 )^2 / γ
γ = norm( s0 )^2

p0 = s0 + β * p0

ai


 