#################################################################
#                   Includes and Modules                        #
#################################################################

using LinearAlgebra, Arpack, Random, SparseArrays
using BenchmarkTools, Printf, TimerOutputs
Random.seed!(12345)

include("../m1Accel.jl")

include("../HOLaGraf_lsmr.jl")
using .HOLaGraf_lsmr
include("../generateDelauney.jl")

# Plotting thingies 
using Plots, ColorSchemes, Plots.PlotMeasures,  LaTeXStrings
pgfplotsx()
theme(:mute)
Plots.scalefontsizes(1.25)
cols=ColorSchemes.Spectral_11;


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
%TODO: REDO !!!
"""
function checkOrder(order, edg2Trig, Δ)
      m = size( order, 1 );

      usedTrigs = Set(1 : Δ);
      flags = falses( m - 1, 1 );
      final = 0;
      for i in 1 : m - 1
            flags[ i ] = length(
                  intersect( usedTrigs, edg2Trig[ order[ i ] ], edg2Trig[ order[ i + 1 ] ] )
            ) > 0;
            setdiff!( usedTrigs, edg2Trig[ order[ i ] ] );
            if isempty(usedTrigs)
                  final = i;
            end
      end
      return flags, sum( flags ), final
end

"""
    randomOrderChecker(B2; maxTries = 5000)

      #TODO:
"""
function randomOrderChecker(edg2Trig, Δ; maxTries = 5000)
      m = size(edg2Trig, 1);
      for repeat in 1 : maxTries
            order = shuffle(1:m);
            flags, s, final = checkOrder(order, edg2Trig, Δ);
            if ( s == 0 ) || ( ( s == 1 ) && ( flags[final] ) )
                  return repeat, order
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



#################################################################
#            Functions of Ordering Timing                       #
#################################################################

"""
    orderGraph(  )

      #TODO:
"""
function getTimedGraphOrdering( G::NiceGraph, edg2Trig, trig2Edg )
      times = @elapsed begin
            trigSet = Set(1:size(G.trigs, 1));
            edgeSet = Set(1:size(G.edges, 1));
      
            order = [ ];
            previous = 0;
            while true
                  flag, newEdge, trigSet = orderingIteration( trigSet, G.edges, edg2Trig, trig2Edg, previous);
                  
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
      end
            
      Δ = size(G.trigs, 1);
      flags, s, final = checkOrder( order, edg2Trig, Δ )
      if ( s > 0 ) && !( ( s==1 ) && ( flags[final] ) ) 
            return false, order, times,  flags
      end
      return true, order, times, flags 
end


"""
    generateComplConGraph( n::Int64 )

      #TODO:
"""
function generateComplConGraph( n::Int64, rem::Int64 )
      edges, trigs = generateComplete( n );
      w_Δ = rand( size(trigs, 1), 1 ) * 0.75 .+ 0.25;
      edges2, trigs2, w_Δ2 = massRemove( rem, edges, trigs, w_Δ, n );

      edg2Trig = getEdges2Trig( edges2, trigs2 );
      trig2Edg = getTrig2Edg( edges2, trigs2, edg2Trig );

      w2 = ones(size(edges2, 1), 1);
      G = NiceGraph( n, edges2, trigs2, w2, w_Δ2 );
      return G, edg2Trig, trig2Edg
end


#################################################################
#            End of Functions of Ordering Timing                #
#################################################################




# C_n - n random edges:


rep = 20;
Ns = [7; 10; 15; 20; 30; 40; 50];
times = zeros( rep, size(Ns, 1) );
maxTrigAdj = zeros( rep, size(Ns, 1) );
aveTrigAdj = zeros( rep, size(Ns, 1) );
trigNum = zeros( rep, size(Ns, 1) );

for i in eachindex(Ns)
      for j in 1:rep
            G, edg2Trig, trig2Edg = generateComplConGraph( Ns[ i ], Ns[ i ] );
            maxTrigAdj[ j, i ] = maximum( length.(edg2Trig) );
            aveTrigAdj[ j, i ] = mean( length.(edg2Trig) );
            trigNum[ j, i ] = size( G.trigs, 1 );
            check, order, times[ j, i ], flags = getTimedGraphOrdering( G, edg2Trig, trig2Edg );
      end
end

# completely connected graph

times2 = zeros( rep, size(Ns, 1) );
maxTrigAdj2 = zeros( rep, size(Ns, 1) );
aveTrigAdj2 = zeros( rep, size(Ns, 1) );
trigNum2 = zeros( rep, size(Ns, 1) );

for i in eachindex(Ns)
      for j in 1:rep
            G, edg2Trig, trig2Edg = generateComplConGraph( Ns[ i ], 0 );
            maxTrigAdj2[ j, i ] = maximum( length.(edg2Trig) );
            aveTrigAdj2[ j, i ] = mean( length.(edg2Trig) );
            trigNum2[ j, i ] = size( G.trigs, 1 );
            check, order, times2[ j, i ], flags = getTimedGraphOrdering( G, edg2Trig, trig2Edg );
      end
end


# triangulation

times3 = zeros( rep, size(Ns, 1) );
maxTrigAdj3 = zeros( rep, size(Ns, 1) );
aveTrigAdj3 = zeros( rep, size(Ns, 1) );
trigNum3 = zeros( rep, size(Ns, 1) );

for i in eachindex(Ns)
      for j in 1:rep
            points, edges2, trigs2 = generateDelauney( Ns[ i ] );
            w_Δ = rand( size(trigs2, 1), 1 ) * 0.75 .+ 0.25;
            edg2Trig = getEdges2Trig( edges2, trigs2 );
            trig2Edg = getTrig2Edg( edges2, trigs2, edg2Trig );
            w2 = ones(size(edges2, 1), 1);
            G = NiceGraph( Ns[ i ] + 4, edges2, trigs2, w2, w_Δ );

            maxTrigAdj3[ j, i ] = maximum( length.(edg2Trig) );
            aveTrigAdj3[ j, i ] = mean( length.(edg2Trig) );
            trigNum3[ j, i ] = size( G.trigs, 1 );
            check, order, times3[ j, i ], flags = getTimedGraphOrdering( G, edg2Trig, trig2Edg );
      end
end




#################################################################
#                     Plotting of rdering Timing                #
#################################################################

# poly
using Polynomials


l = @layout [ grid( 1, 2 ) ]

plot( layout = l,  legend=:topleft )
res = Polynomials.fit( log10.(Ns), vec( log10.(median(times; dims = 1)') ), 1 );
plot!( Ns, median(times; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4), sp = 1,
      color = cols[ 9 ];
      marker = (:diamond, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"C_n - n \textrm{ \; random edges}, \gamma = %$( round( res.coeffs[2] ; digits = 2) ) "
)
res = Polynomials.fit( vec(log10.(median(trigNum; dims = 1)' .* median(aveTrigAdj; dims = 1)')), vec( log10.(median(times; dims = 1)') ), 1 );
plot!( median(trigNum; dims = 1)' .* median(aveTrigAdj; dims = 1)', median(times; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4), sp = 2,
      color = cols[ 9 ];
      marker = (:diamond, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"C_n - n \textrm{ \; random edges}, \gamma = %$( round( res.coeffs[2] ; digits = 2) ) "
)
res = Polynomials.fit( log10.(Ns), vec( log10.(median(times2; dims = 1)') ), 1 );
plot!( Ns, median(times2; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4),  sp = 1,
      color = cols[ 1 ];
      marker = (:circle, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"C_n  \textrm{ \; random edges}, \gamma = %$( round( res.coeffs[2] ; digits = 2) )"
)
res = Polynomials.fit( vec(log10.(median(trigNum2; dims = 1)' .* median(aveTrigAdj2; dims = 1)')), vec( log10.(median(times2; dims = 1)') ), 1 );
plot!( median(trigNum2; dims = 1)' .* median(aveTrigAdj2; dims = 1)', median(times2; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4),  sp = 2,
      color = cols[ 1 ];
      marker = (:circle, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"C_n  \textrm{ \; random edges}, \gamma = %$( round( res.coeffs[2] ; digits = 2) ) "
)
res = Polynomials.fit( log10.(Ns), vec( log10.(median(times3; dims = 1)') ), 1 );
plot!( Ns .+ 4, median(times3; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4),  sp = 1,
      color = cols[ 4 ];
      marker = (:square, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"\textrm{ triangulation}, \gamma = %$( round( res.coeffs[2] ; digits = 2) )"
)
res = Polynomials.fit( vec(log10.(median(trigNum3; dims = 1)' .* median(aveTrigAdj3; dims = 1)')), vec( log10.(median(times3; dims = 1)') ), 1 );
plot!( median(trigNum3; dims = 1)' .* median(aveTrigAdj3; dims = 1)', median(times3; dims = 1)', 
      xscale = :log10, yscale = :log10,  
      line = (:dot, 4),  sp = 2,
      color = cols[ 4 ];
      marker = (:square, 6, 0.8, Plots.stroke(1, :black)),
      labels = L"\textrm{ triangulation}, \gamma = %$( round( res.coeffs[2] ; digits = 2) ) "
)
xticks!( Ns, latexstring.(Ns) , sp = 1)
xlabel!(L"| \mathcal V_0 ( \mathcal K ) |",  sp = 1)
xlabel!(L"| \mathcal V_2 ( \mathcal K ) | \cdot \overline{\Delta}_e",  sp = 2)
ylabel!(L"\textrm{exec \; time}",  sp = 1)
title!(L"\textrm{Ordering \; time \; vs. }n",  sp = 1)
title!(L"\textrm{Ordering \; time \; vs. triangles}",  sp = 2)
plot!( size = (1000, 400) )


savefig("order_complexity/order.tex")
savefig("order_complexity/order.pdf")







