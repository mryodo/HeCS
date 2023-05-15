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


constructEdgeMeasure( w, edge2Trig ) = [ sum( w[collect(edge2Trig[i]) ] ) for i in 1 : size(edge2Trig, 1) ] 


n, edges, trigs = generateComplete(6)
w = 0.75*rand( size(trigs, 1), 1 ) .+ 0.25

edge2Trig = getEdges2Trig( edges, trigs )
trig2Edge = getTrig2Edg( edges, trigs, edge2Trig )

p = constructEdgeMeasure( w, edge2Trig )
σ = wsample( Array(1 : size(edge2Trig, 1)), p )






function intersectTriangles( t1, t2, trig2Edge, edge2Trig )
      common = intersect( trig2Edge[t1], trig2Edge[t2] )
      if !isempty(common)
            σ = collect(common)[1]
            if !isempty( edge2Trig[ σ ] )
                  return true
            end
      end
      return false
end


extractAdjacency( t, w, edge2Trig, trig2Edge ) =  sum( [ sum( [ 1/( 1 / w[t2]+ 1 / w[t] ) for t2 in edge2Trig[σ] ] )  for σ in trig2Edge[t]  ] ) 


function constructTrigMeasure( σ, w, edge2Trig, trig2Edge ) 
      q = zeros( length(edge2Trig[ σ ]) ) 
      val = zeros( length(edge2Trig[ σ ]) ) 
      tmp = collect(edge2Trig[ σ ])
      for i in 1:length(edge2Trig[ σ ] )
            q[i]= extractAdjacency( tmp[i] , w, edge2Trig, trig2Edge )  
            val[i] = tmp[i]
      end
      return q, Int.(val)
end


p = constructEdgeMeasure( w, edge2Trig )
σ = wsample( Array(1 : size(edge2Trig, 1)), p )
q, val = constructTrigMeasure( σ, w, edge2Trig, trig2Edge )
τ = wsample( val, q )

B2 = B2fromTrig(edges, trigs)
C = Array{Float64}(undef, size(edges, 1), 0)

p = constructEdgeMeasure( w, edge2Trig )
σ = wsample( Array(1 : size(edge2Trig, 1)), p )
q, val = constructTrigMeasure( σ, w, edge2Trig, trig2Edge )
τ = wsample( val, q )

C = [C sqrt(w[τ])*B2[:, τ] ]

edge2Trig_new = [ setdiff(edge2Trig[i], edge2Trig[σ])
      for i in 1:length(edge2Trig)
]


p = constructEdgeMeasure( w, edge2Trig_new )
σ = wsample( Array(1 : size(edge2Trig_new, 1)), p )
q, val = constructTrigMeasure( σ, w, edge2Trig_new, trig2Edge )
τ = wsample( val, q )

C = [C sqrt(w[τ])*B2[:, τ] ]
edge2Trig_new = [ setdiff(edge2Trig_new[i], edge2Trig_new[σ])
      for i in 1:length(edge2Trig_new)
]


τ=10
σ = 9
C = [C sqrt(w[τ])*B2[:, τ] ]
edge2Trig_new = [ setdiff(edge2Trig_new[i], edge2Trig_new[σ])
      for i in 1:length(edge2Trig_new)
]


p = constructEdgeMeasure( w, edge2Trig_new )
σ = wsample( Array(1 : size(edge2Trig_new, 1)), p )
q, val = constructTrigMeasure( σ, w, edge2Trig_new, trig2Edge )
τ = wsample( val, q )

C = [C sqrt(w[τ])*B2[:, τ] ]
edge2Trig_new = [ setdiff(edge2Trig_new[i], edge2Trig_new[σ])
      for i in 1:length(edge2Trig_new)
]


τ=2
σ = 1
C = [C sqrt(w[τ])*B2[:, τ] ]
edge2Trig_new = [ setdiff(edge2Trig_new[i], edge2Trig_new[σ])
      for i in 1:length(edge2Trig_new)
]


τ=7
σ = 5
C = [C sqrt(w[τ])*B2[:, τ] ]
edge2Trig_new = [ setdiff(edge2Trig_new[i], edge2Trig_new[σ])
      for i in 1:length(edge2Trig_new)
]


LU = B2 * diagm(vec(w)) * B2'
function condPlus( A; thr = 1e-8 )
      m = size(A, 1);
      #σ = svds(L1up, nsv = m - 1)[1].S;
      σ = svd(Matrix(A)).S;
      return maximum( σ ) / minimum( σ[ abs.(σ) .> thr ])
end

C = Matrix(C)

condPlus(LU)
condPlus(pinv(C)* Π' * LU * Π * pinv(C'))
condPlus(pinv( Π * C) * LU * pinv(C' * Π' ))

perm = [7; 4; 9; 2; 1; 5; 3; 6; 8; 10]
Π = zeros( size( edges, 1 ), size(edges, 1) )

for i in 1:size(edges, 1)
      Π[ i, perm[i] ] = 1
end
 

function getTrigAdj( edges, trigs, w, edge2Trig, trig2Edge )
      Δ = size(trigs, 1)
      A = zeros( Δ, Δ )

      for i in 1 : Δ-1
            for j in i+1 : Δ
                  if !isempty( intersect( trig2Edge[i], trig2Edge[j] ) )
                        A[ i, j ] = 1 / ( 1 / w[i] + 1 / w[j] )
                        A[ j, i ] = 1 / ( 1 / w[i] + 1 / w[j] )
                  end
            end
      end
      return A
end

A = getTrigAdj( edges, trigs, w, edge2Trig, trig2Edge )
deg = A * ones( size(trigs, 1) )

W = diagm( vec(w) )
res = diag( sqrt.(W) * B2' * pinv( B2 * W * B2' ) * B2 * sqrt.(W) )

Δ = size( trigs, 1 )
num = Int( round(2/3 * size(edges, 1) )) - 1

new = wsample( Array(1 : Δ), deg, num; replace = false )




















