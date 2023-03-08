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

N = 10;   # number of additional points in the grid
n  = N + 4;

points, edges, trigs = generateDelauney( N )

ν_init = size(edges, 1 ) / ( ( N + 4) * ( N + 3 ) / 2 );
ν = 0.3;
#ν = ν_init;

backlash = Int( - size( edges, 1 ) + round( ν * n * (n-1) / 2 ) ) 

if backlash < 0
      for i in 1 : -backlash
            global edges, n, trigs
            ind = getIndx2Kill( edges ) ;
            edges, trigs = killEdge(ind, n, edges, trigs)
      end
else
      for i in 1 : backlash
            global edges, n, trigs
            ind = getNewEdge(n, edges);
            edges, trigs = addEdge(ind, n, edges, trigs)
      end
end


edge2Trig = getEdges2Trig(edges, trigs)
trig2Edg = getTrig2Edg(edges, trigs, edge2Trig)

Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

Σ = [];

while !isempty(Free)
      global Σ, edge2Trig, trig2Edg, Free, Ls
      σ = pop!(Free) 
      Σ = [ Σ; σ ]
      τ = first( edge2Trig[ σ ] ) 
      for ξ in trig2Edg[ τ ]
            setdiff!( edge2Trig[ ξ ], Set([ τ ]) )
            Ls[ ξ ] = Ls[ ξ ] - 1
            ( Ls[ ξ ] == 1 ) && ( union!(Free, Set( [ ξ ] ) ) )
            ( Ls[ ξ ] == 0 ) && ( setdiff!(Free, Set( [ ξ ] ) ) )
      end
end

sum( Ls )