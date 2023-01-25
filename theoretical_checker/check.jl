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
function intersectEdges(edge1, edge2, edg2Trig)
      #return length( intersect( Set(edge1), Set(edge2) ) ) > 0 # OLD  VERSION: do not work
      return length( intersect( edg2Trig[edge1], edg2Trig[edge2] ) ) > 0
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
                  #TODO: FIX  THIS SHIT NOW
                  if ( previous == 0 ) || !intersectEdges( previous, newEdge, edg2Trig)  #!intersectEdges( edges2[previous, :], edges2[newEdge, :] )
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


"""
    getStar( trigList, G )

TBW
"""
function getStar( trigList, G )
      m = size( G.edges, 1 );
      A = zeros( m, m );
      for t in trigList
            A = A + G.w_Δ[t] * G.B2[ : , t] * G.B2[ : , t]';
      end
      return A
end

function getCyclic( trigList, G )
      m = size( G.edges, 1 );
      K = zeros( m, m );
      Ω = sum( [ G.w_Δ[t] for t in trigList ] )
      for t1 in trigList
            for t2 in trigList
                  K = K + G.w_Δ[t1] * G.w_Δ[t2] * ( G.B2[:, t1] - G.B2[:, t2]) * ( G.B2[:, t1] - G.B2[:, t2])'; 
            end
      end
      return 1 / (2 * Ω) * K
end

function δ(i)
      δ = zeros(m, 1); δ[i] = 1.0;
      return δ
end


function Schur(Si, i)
      return Si - 1 / ( δ(i)' * Si * δ(i) )[1, 1] * Si * δ(i) * δ(i)' * Si 
end



N = 5;
G, edg2Trig, trig2Edg = generateComplConGraph( N, 0 );
m = size( G.edges, 1 );
Δ = size( G.trigs, 1 );

check, order, times, flags = getTimedGraphOrdering( G, edg2Trig, trig2Edg );
G, edg2Trig, trig2Edg = reodering( order, G, edg2Trig );
check

ε = 1e-8
G.w_Δ[ 3 ] = ε; G.w_Δ[ 7 ] = ε;



L1up = getL1up(G);

vals = eigen(Matrix(L1up)).values




usedTrigs = Set( 1:Δ )
S0 = getStar( usedTrigs, G )

usedTrigs1 = setdiff( usedTrigs, edg2Trig[1] )
S1 = getStar( usedTrigs1, G )
K1 = getCyclic( edg2Trig[1], G)

usedTrigs2 = setdiff( usedTrigs1, edg2Trig[2] )
S2 = getStar( usedTrigs2, G )
K2 = getCyclic( setdiff(edg2Trig[2], edg2Trig[1]), G)

usedTrigs3 = setdiff( usedTrigs2, edg2Trig[3] )
S3 = getStar( usedTrigs3, G )
K3 = getCyclic( setdiff( setdiff(edg2Trig[3], edg2Trig[1]), edg2Trig[2] ) , G)

usedTrigs4 = setdiff( usedTrigs3, edg2Trig[4] )
S4 = getStar( usedTrigs4, G )
K4 = getCyclic( setdiff( setdiff( setdiff(edg2Trig[4], edg2Trig[1]), edg2Trig[2] ), edg2Trig[3] ) , G)

usedTrigs5 = setdiff( usedTrigs4, edg2Trig[5] )
S5 = getStar( usedTrigs5, G )



E0 = L1up

i = 1
E1 = Schur(E0, i)

i = 2
E2 = Schur(E1, i)

i = 3
E3 = Schur(E2, i)

i = 4
E4 = Schur(E3, i)

i = 5
E5 = Schur(E4, i)





# restrictive set
σ_i = Set([]);

Si = L1up;
usedTrig = Set( 1 : Δ );
Si_log = zeros( m, m, m );
C = zeros( m , m) ;
for i in 1 : m 
      global Si
      if abs(Si[ i , i ]) < 1e-8 
            break
      end
      C[:, i ] = 1 / sqrt(Si[ i , i ]) * Si[ : , i ];
      Si = Si - 1 / Si[ i , i ] * Si[ : , i ] * Si[ :, i]';
      Si_log[ :, :, i] = Si;    
end
sum( abs.(L1up - C * C') .> 1e-8 ) 



C2 = zeros( m , m) ;
σ_i = Set([]) ;
usedTrig = Set( 1 : Δ ) ;
for i in 1:m
      global σ_i, usedTrig, G, edg2Trig, C2
      (isempty(usedTrig)) && break;
      local A = zeros( m, m );
      for t in usedTrig
            A = A + G.w_Δ[ t ] * G.B2[ : , t ] * G.B2[ : , t]';
      end
      C2[:, i ] = 1 / sqrt(A[ i , i ]) * A[ : , i ];
      setdiff!( usedTrig, edg2Trig[ i ] ); 
end
sum( abs.(L1up - C2 * C2') .> 1e-8 ) 
eigs( pinv(C2) * L1up * pinv(C2)')[1]


usedTrig = Set( 1 : Δ ) ;
A = zeros( m, m );
for t in usedTrig
      A = A + G.w_Δ[ t ] * G.B2[ : , t ] * G.B2[ : , t]';
end


for i in 1 : m-1
      global Si, usedTrig, σ_i, G, edg2Trig
      # exact Schur decomp S_i 
      exactSi = Si - 1 / Si[ i , i  ] * Si[ : , i ] * Si[ : , i ]'; 
      exactCi = 1 / sqrt( exactSi[i , i] ) * exactSi[ : , i ];

      union!(σ_i, Set([i]));
      setdiff!(usedTrig, edg2Trig[i] );

      # star-term of Si
      A = zeros( m, m );
      for t in usedTrig
            A = A + G.w_Δ[t] * G.B2[:, t ] * G.B2[: , t ]' ;
      end

      starSi = A;
      starCi = 1 / sqrt( starSi[i , i] ) * starSi[:, i];
       
      # clique term
      tmp =  [ G.w_Δ[ t ] for t in edg2Trig[ i ] ]
      isempty(tmp) ? Ω = 0 :  Ω = sum( tmp );
      K = zeros( m, m);
      for t1 in edg2Trig[i]
            for t2 in edg2Trig[i]
                  K = K + G.w_Δ[t1] * G.w_Δ[t2] * ( G.B2[ : , t1 ] - G.B2[ :, t2 ] ) * ( G.B2[:, t1 ] - G.B2[:, t2] )';
            end
      end
      K = 1 / ( 2 * Ω ) * K;

      for j in i+1:m
            setdiff!(edg2Trig[j], edg2Trig[i])
      end

      Si = A + K ;
      Ci = 1 / sqrt( Si[i , i] ) * Si[ :, i];

      #@printf "error: %f \n" norm(Si-exactSi);
      #@printf " c error: %f %f \n\n" norm(exactCi - Ci) norm(starCi-Ci)

end
























