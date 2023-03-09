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

      while !isempty(Free)
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
      return sum( Ls ) == 0, Σ, edge2Trig, trig2Edg, Ls, Free  
end





Random.seed!(0)

n, points, ν_init, edges, trigs = sparseDelaunay( N = 10, ν = 0.4)
flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )




N = 20
ν_init = ( 3*(N+4) - 7 ) / ( (N+4)*(N+3)/2 )
νs = Array( ν_init*1.005:0.001: ν_init*1.25 )
rep = 500
freq = zeros( size(νs, 1), 1 )

for i in 1:size(νs, 1)
      global νs, rep, N, freq
      ν = νs[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, ν_init, edges, trigs = sparseDelaunay( N = N, ν = ν)
            flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
            !flag && ( freq[i] = freq[i] + 1 )
      end

end

freq = freq / rep



N = 10
ν_init01 = ( 3*(N+4) - 7 ) / ( (N+4)*(N+3)/2 )
νs01 = Array( ν_init01*1.005:0.001: ν_init01*1.25 )
rep = 500
freq01 = zeros( size(νs01, 1), 1 )

for i in 1:size(νs01, 1)
      global νs01, rep, N, freq01
      ν = νs01[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, ν_init01, edges, trigs = sparseDelaunay( N = N, ν = ν)
            flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
            !flag && ( freq01[i] = freq01[i] + 1 )
      end

end

freq01 = freq01 / rep

N = 15
ν_init015 = ( 3*(N+4) - 7 ) / ( (N+4)*(N+3)/2 )
νs015 = Array( ν_init015*1.005:0.001: ν_init015*1.25 )
rep = 500
freq015 = zeros( size(νs015, 1), 1 )

for i in 1:size(νs015, 1)
      global νs015, rep, N, freq015
      ν = νs015[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, ν_init015, edges, trigs = sparseDelaunay( N = N, ν = ν)
            flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
            !flag && ( freq015[i] = freq015[i] + 1 )
      end

end

freq015 = freq015 / rep

plot()

plot!( νs[vec(freq.>0)]/ν_init, freq[vec(freq.>0)], marker = :circle, xscale=:log10, yscale=:log10, label=L"\mathcal{V}_0(\mathcal K)=24")
#plot!( [ν_init, ν_init], [0, 1], color=:black )

plot!( νs01[vec(freq01.>0)]/ν_init01, freq01[vec(freq01.>0)], marker = :circle, xscale=:log10, label=L"\mathcal{V}_0(\mathcal K)=14" )
#plot!( [ν_init01, ν_init01], [0, 1], color=:black )

plot!( νs015[vec(freq015.>0)]/ν_init015, freq015[vec(freq015.>0)], marker = :circle, xscale=:log10, label=L"\mathcal{V}_0(\mathcal K)=19")
#plot!( [ν_init015, ν_init015], [0, 1], color=:black )

plot!([1, 1], [0.25, 1], color=:black, linestyle=:dash, label="")
xticks!( [1, 1.1, 1.25], [L"\nu_\Delta", L"1.1\nu_\Delta", L"1.25\nu_\Delta"], tickfontsize=15 )
yticks!( [ 0.25, 0.5, 1], [L"0.25", L"0.5", L"1"])
xlabel!(L"\mathrm{sparsity \; (wrt \; triangulation)}")
ylabel!(L"\mathrm{probability \; of \; 2-Core}")
plot!( size=( 600, 600 ), legend=:bottomright, legendfont=font(18) )

savefig("triang_prob.tex")
savefig("triang_prob.pdf")





#                                                                      #

dist( x, y ) = norm( x - y, 2 )

function generateSensor( N = 10 , eps = 0.5)
      points =  rand(N, 2)
      edges = Array{Integer}(undef, 0, 2)
      for i in 1 : N-1
            for j in i+1 : N
                  if dist( points[i, :], points[j, :] ) < eps
                        edges = [edges; i j ]
                  end
            end
      end
      edges

      trigs = Array{Integer}(undef, 0, 3)

      for i in 1 : size( edges, 1 )-1
            for j in i+1 : size( edges, 1 )
                  if edges[ i, 1 ] == edges[ j, 1 ]
                        if sum( all( edges .== sort([ edges[i, 2] edges[j, 2] ]; dims=2)  , dims = 2) ) > 0
                              trigs = [ trigs; sort( [ edges[i, 1] edges[i, 2] edges[j, 2] ]; dims = 2 ) ]
                        end
                  end
            end
      end
      return N, points, edges, trigs
end

N = 20
eps = 0.5

Random.seed!(1)
n, points, edges, trigs = generateSensor(N, eps)
flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )





N = 20
ϵs = Array(0.1:0.005:0.45)
rep = 500
freq = zeros( size(ϵs, 1), 1 )

for i in 1:size(ϵs, 1)
      global ϵs, rep, N, freq
      eps = ϵs[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, edges, trigs = generateSensor(N, eps)
            if size(edges, 1) > 0
                  flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
                  !flag && ( freq[i] = freq[i] + 1 )
            end
      end

end

freq = freq / rep


N = 10
ϵs01 = Array(0.1:0.005:0.45)
rep = 500
freq01 = zeros( size(ϵs01, 1), 1 )

for i in 1:size(ϵs01, 1)
      global ϵs01, rep, N, freq01
      eps = ϵs01[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, edges, trigs = generateSensor(N, eps)
            if size(edges, 1) > 0
                  flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
                  !flag && ( freq01[i] = freq01[i] + 1 )
            end
      end

end

freq01 = freq01 / rep



N = 15
ϵs015 = Array(0.1:0.005:0.45)
rep = 500
freq015 = zeros( size(ϵs015, 1), 1 )

for i in 1:size(ϵs015, 1)
      global ϵs015, rep, N, freq015
      eps = ϵs015[i]
      for j in 1:rep
            Random.seed!(j)
            n, points, edges, trigs = generateSensor(N, eps)
            if size(edges, 1) > 0
                  flag, Σ, edge2Trig, trig2Edg, Ls, Free = greedyCollapse( edges, trigs )
                  !flag && ( freq015[i] = freq015[i] + 1 )
            end
      end
end

freq015 = freq015 / rep




plot()
plot!( ϵs, freq, marker = :circle, label = L"\mathcal{V}_0(\mathcal K) = 20" )
plot!( ϵs01, freq01, marker = :circle, label = L"\mathcal{V}_0(\mathcal K) = 10" )
plot!( ϵs015, freq015, marker = :circle, label = L"\mathcal{V}_0(\mathcal K) = 15" )
xlabel!(L"\mathrm{percolation, \;} \varepsilon")
ylabel!(L"\mathrm{probability \; of \; 2-Core}")
plot!(size=(600, 600), legend = :bottomright, legendfont = font(18) )

savefig("percol_prob.tex")
savefig("percol_prob.pdf")
















N = 20

eps = 0.2

edges = Array{Integer}(undef, 0, 2)
for i in 1 : N-1
      for j in i+1 : N
            if dist( points[i, :], points[j, :] ) < eps
                  edges = [edges; i j ]
            end
      end
end
edges

trigs = Array{Integer}(undef, 0, 3)

for i in 1 : size( edges, 1 )-1
      for j in i+1 : size( edges, 1 )
            if edges[ i, 1 ] == edges[ j, 1 ]
                  if sum( all( edges .== sort([ edges[i, 2] edges[j, 2] ]; dims=2)  , dims = 2) ) > 0
                        trigs = [ trigs; sort( [ edges[i, 1] edges[i, 2] edges[j, 2] ]; dims = 2 ) ]
                  end
            end
      end
end













