#################################################################
#                   Includes and Modules                        #
#################################################################

using LinearAlgebra, Arpack, Random, SparseArrays
using BenchmarkTools, Printf, TimerOutputs
Random.seed!(12345)

include("m1Accel.jl")

include("HOLaGraf_lsmr.jl")
using .HOLaGraf_lsmr
include("generateDelauney.jl")

#= Plotting thingies: MOVE TO https://github.com/liuyxpp/MakiePublication.jl.git #TODO:

using Plots, ColorSchemes, Plots.PlotMeasures,  LaTeXStrings
pgfplotsx()
theme(:mute)
Plots.scalefontsizes(1.25)
cols=ColorSchemes.Spectral_11;
=#

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
#                   Functions of Order                          #
#################################################################

"""
    intersectEdges(edge1, edge2)

check if two edges intersect
"""
function intersectEdges(edge1, edge2)
      return length( intersect( Set(edge1), Set(edge2) ) ) > 0
end

"""
    orderingIteration( trigSet, edges2, trig2Edg, previous = 0 )

      #TODO:
"""
function orderingIteration( trigSet, edges2, edg2Trig, trig2Edg, previous = 0 )
      tmpTrigSet = deepcopy( trigSet );
      
      while !isempty(tmpTrigSet)
            newTrig = first( tmpTrigSet );
            edgList = collect(trig2Edg[newTrig]);
            i = 1;
            newEdge = 0;
            while i < 4
                  newEdge = edgList[i]; 
                  if ( previous == 0 ) || !intersectEdges( edges2[previous, :], edges2[newEdge, :] )
                        break
                  end
                  i = i + 1;
            end
            if i < 4
                  setdiff!( trigSet, edg2Trig[newEdge] )
                  return true, newEdge, trigSet
            else
                  delete!(tmpTrigSet, newTrig);
            end
      end
      return false, 0, trigSet
end

"""
    checkOrder(order, B2)

For the matrix `B2` and the vector `order` checks for every 2 consequetive triangles in `order` if there and intersection.
"""
function checkOrder(order, edg2Trig)
      m = size( order, 1 );
      flags = [ length( intersect( edg2Trig[order[i]], edg2Trig[order[i+1]] ) ) for i in 1 : m - 1 ];
      return flags, sum( flags )
end

"""
    randomOrderChecker(B2; maxTries = 5000)

      #TODO:
"""
function randomOrderChecker(B2; maxTries = 5000)
      m = size(B2, 1);
      for repeat in 1 : maxTries
            order = shuffle(1:m);
            flags, s = checkOrder(order, B2);
            if ( s == 0 ) || ( ( s == 1 ) && ( flags[end] = 1 ) )
                  return order
            end
      end
end

"""
    reodering( order, G, edg2Trig, trig2Edg )

      #TODO:
"""
function reodering( order, G, edg2Trig )

      G.edges = G.edges[order, :];
      G.B1 = G.B1[:, order];
      G.B2 = G.B2[order, :];
      G.w = G.w[order];
      edg2Trig = edg2Trig[ order ];
      trig2Edg = getTrig2Edg( G.edges, G.trigs, edg2Trig );

      for i in axes(G.B2, 2)
            ind = findfirst(.!(G.B2[:, i] .== 0));
            if G.B2[ind, i] < 0
                  G.B2[:, i] = (-1.0) * G.B2[:, i];
            end
      end

      return G, edg2Trig, trig2Edg
end


#################################################################
#               END of Functions of Order                       #
#################################################################


#################################################################
#                   Functions of Graph Generation               #
#################################################################

"""
    generateComplete( n )

      #TODO:
"""
function generateComplete( n )
      edges = zeros(Int, Int(n*(n-1)/2), 2);
      trigs = zeros(Int, Int(n*(n-1)*(n-2)/6), 3);

      cur = 0; 
      for i in 1:n-1
            for j in i + 1:n
                  cur = cur + 1;
                  edges[ cur , :] = [i j];
            end
      end

      cur = 0; 
      for i in 1:n-2
            for j in i + 1:n-1
                  for k in j + 1:n
                        cur = cur + 1;
                        trigs[ cur , :] = [i j k];
                  end
            end
      end
      return edges, trigs
end
 
"""
    massRemove( remove, edges, trigs, w_Δ, n )

      #TODO:
"""
function massRemove( remove, edges, trigs, w_Δ, n )
      edges2 = edges;
      trigs2 = trigs;
      w_Δ2 = w_Δ;
      
      for i in 1 : remove
            indx = getIndx2Kill(edges2);
            edges2, trigs2, w_Δ2 = killEdge(indx, n , edges2, trigs2, w_Δ2);
      end

      return edges2, trigs2, w_Δ2
end

#################################################################
#            End of Functions of Graph Generation               #
#################################################################

n = 6;
edges, trigs = generateComplete( n );
w_Δ = rand( size(trigs, 1), 1 ) * 0.75 .+ 0.25;
edges2, trigs2, w_Δ2 = massRemove( n, edges, trigs, w_Δ, n );

edg2Trig = getEdges2Trig( edges2, trigs2 );
trig2Edg = getTrig2Edg( edges2, trigs2, edg2Trig );

w2 = ones(size(edges2, 1), 1);
G = NiceGraph( n, edges2, trigs2, w2, w_Δ2 );

trigSet = Set(1:size(trigs2, 1));
edgeSet = Set(1:size(edges2, 1));

order = [ ];
previous = 0;
while true
      global trigSet, edgeSet, G, edges2, trigs2, edg2Trig, trig2Edg, order, previous

      flag, newEdge, trigSet = orderingIteration( trigSet, edges2, edg2Trig, trig2Edg, previous);
      
      if !flag
            break
      end
      order = [ order; newEdge ];
      previous = newEdge;
      delete!(edgeSet, previous);
end

if !isempty(edgeSet)
      order = [ order; collect( edgeSet ) ];
end

checkOrder( order, edg2Trig )






















