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
using Distributions
import Base.Iterators: flatten
using Polynomials
using Krylov
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

function greedyCollapseShort( edge2Trig, trig2Edg, perm, edges, iperm, Ls, Free)
      
      #Ls = [ length(edge2Trig[i]) for i in 1 : size( edges, 1 ) ]
      #Free = Set([ i for i in 1 : size(edges, 1 ) if Ls[i] == 1 ])

      Σ = [];
      Τ = [ ]; 
      while !isempty(Free)
            σ = pop!(Free) 
            Σ = [ Σ; σ ]
            τ = first( edge2Trig[ σ ] ) 
            Τ = [ Τ ; τ ]
            #in = findall(perm .== τ)[1][1] 
            #in = iperm[τ]
            in = τ
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

function B1fromEdges(n, edges)
      m = size(edges, 1);
      B1 = spzeros(n, m);
      
      for i in 1:m
          B1[edges[i, 1], i] = -1;
          B1[edges[i, 2], i] = 1;
      end
      return B1
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



function cgls( A, b; tol = 1e-3, maxit = 500 )
      x0 = zeros( size(A, 2), 1 )
      
      r0 = b - A * x0
      p0 = A' * r0
      s0 = p0

      γ = norm(s0)^2
      β = 2.0
      it = 0 
      for i in 1 : maxit
            it = i
            qi = A * p0
            ai = γ / norm(qi)^2
            x0 = x0 + ai * p0
            r0 = r0 - ai * qi

            if norm(r0, Inf) < tol
                  break
            end

            s0 = A' * r0
            β = norm( s0 )^2 / γ
            γ = norm( s0 )^2

            p0 = s0 + β * p0
      end

      return x0, it 
end

function getAllEdges( n )
      allEdges = Array{Integer}(undef, 0, 2)
      for i in 1:(n-1) 
            for j in (i+1):n
                  allEdges = [ allEdges; i j  ]
            end
      end
      return allEdges
end

function repeatTries(N, add, rep)
      κ_original = zeros( rep )
      κ_precon = zeros( rep )
      κ_chol = zeros( rep )
      κ_ldl = zeros( rep )
      it_ldl = zeros( rep )
      it_original = zeros( rep )
      it_precon = zeros( rep )
      m_sizes = zeros( rep )
      Δ_sizes = zeros( rep )
      times = zeros( rep )
      times_chol = zeros( rep )
      times_precon = zeros( rep )

      for repIt in 1 : rep
            points, edges, trigs = generateDelauney( N )
            edges2, trigs2 = deepcopy(edges), deepcopy(trigs)
            n = N + 4

            ν_Δ = size(edges, 1) / binomial( n, 2 )

            #add = 60 # how many edges do we need to add
            
            allEdges = getAllEdges(n)
            for i in axes(edges2, 1)
                  indx = findall(all(allEdges .== edges2[i, :]', dims=2))[1][1]
                  allEdges = allEdges[ 1:size(allEdges, 1) .!= indx, : ]
            end
            for i in 1 : add      
                  ind, allEdges = getNewEdge2(n, edges2, allEdges);
                  edges2, trigs2 = addEdge(ind, n, edges2, trigs2)
            end

            B1 = B1fromEdges(n, edges2)

            w = zeros( size(trigs2, 1), 1 )
            m = size( edges2, 1 )
            Δ = size( trigs2, 1 )

            m_sizes[ repIt ] = m
            Δ_sizes[ repIt ] = Δ

            if size(Δ, 1) > 0.95*m/4*log(4*m)
                  break
            end

            edge2Trig = getEdges2Trig( edges2, trigs2 )
            trig2Edge = getTrig2Edg( edges2, trigs2, edge2Trig )

            #w_e = abs.( rand( Cauchy(0, 1) , size(edges2, 1) ) )
            #w_e[ w_e .> 3.5 ] .= 3.5
            w_e = abs.( rand( size(edges2, 1) ) )
            for i in 1 : Δ
                  w[i] = minimum( w_e[ collect( trig2Edge[i] ) ]  )
            end

            #μ1 = 1.0; s1 = 1/3; μ2 = 0.2;  s2 = 0.01;
            #=μ1 = 1.0; s1 = 1/3; μ2 = 0.5;  s2 = 0.1;
            for i in 1 : Δ
                  #global w, trigs2, trigs
                  if sum( all( trigs .== trigs2[i, :]'; dims = 2)) == 1   
                        if rand() > 3/Δ
                              w[ i ] = abs( randn()*s1 + μ1 )
                        else
                              w[ i ] = abs( randn()*s2 + μ2/2 )
                        end 
                  else
                        w[ i ] = abs( randn()*s2 + μ2 )
                  end
            end=#

            W = diagm(vec(sqrt.(w)))
            B2 = B2fromTrig( edges2, trigs2 )
            Lu =  B2 * W * W * B2'

            original = condPlus( Lu )

            perm = sortperm( w, dims =1, rev = true )
            iperm = invperm(vec(perm))

            Ls = [ length(edge2Trig[i]) for i in 1 : size( edges2, 1 ) ]
            Free = Set([ i for i in 1 : size(edges2, 1 ) if Ls[i] == 1 ])
      

            timePeriod = @elapsed begin
                  fl, Σ1, Τ1, edge2Trig, trig2Edge, Ls, Free = greedyCollapseShort( edge2Trig, trig2Edge, vec(1:Δ), edges2, vec(1:Δ), Ls, Free )
            end
            #timePeriod = 0
                  #leftovers = Set(unique(collect(flatten(collect.(edge2Trig)))))
                  leftovers = unique(collect(flatten(collect.(edge2Trig))))
                  leftoversOP = sort(iperm[leftovers])
                  sub = [ ] 
                  ind = 1

            #timePeriod2 = @elapsed begin
                  subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
                  #subTrig2Edg = [ ]
                  timePeriod2 = 0
                  iperm2 = zeros(Int, Δ)
                  for i in eachindex( leftoversOP)
                        ind = leftoversOP[i]
                        iperm2[perm[ind]] = i
                        #subEdg2Trigs = [ Set([ ]) for i in 1 : m ]
                        sub = [ sub; perm[ind] ]
                        #subTrig2Edg = [  trig2Edge[ sub_ind ] for sub_ind in sub ]
                        #subTrig2Edg = [ subTrig2Edg; trig2Edge[ perm[ind] ] ]
                        for sub_ind in sub
                              for e in trig2Edge[ sub_ind ]
                                    union!( subEdg2Trigs[e], sub_ind )
                              end
                        end
                        #for e in trig2Edge[ perm[ind] ]
                        #      union!( subEdg2Trigs[e], perm[ind] )
                        #end

                        Ls = [ length(subEdg2Trigs[i]) for i in 1 : size( edges2, 1 ) ]
                        Free = Set([ i for i in 1 : size(edges2, 1 ) if Ls[i] == 1 ])
                  

                        timePeriod2 += @elapsed begin
                              fl, _, _, _, _, Ls, Free = greedyCollapseShort( subEdg2Trigs, trig2Edge, sub, edges2, iperm2, Ls, Free )
                        end
                              subEdg2Trigs[ Ls .>= 0  ] .= empty.(subEdg2Trigs[ Ls .>= 0  ] )
                              if !fl
                                    pop!(sub)
                                    #pop!(subTrig2Edg)
                                    for e in trig2Edge[ perm[ind] ]
                                          setdiff!( subEdg2Trigs[e], perm[ind] )
                                    end
                              end
                        
                  end
                        

                  filter = [Τ1; sub]
            #end
                  Π = indicator( sort(filter), Δ )
                  trigs3 = trigs2[ sort(filter), :]
                  edge2Trig2 = getEdges2Trig( edges2, trigs3 )
                  trig2Edge2 = getTrig2Edg( edges2, trigs3, edge2Trig2 )

                  fl, Σ, Τ, edge2Trig2, trig2Edg2, Ls, Free = greedyCollapse( edges2, trigs3 )
            
            Σ_full = [ Σ; sort(collect(setdiff( Set(1:m), Set(Σ)))) ]
            Τ_full = [ filter[Τ]; sort(collect(setdiff( Set(1:Δ), Set(filter[Τ]))))  ]

            P1 = Perm( Σ_full )
            P2 = Perm( Τ_full )
          
            C = P1 * B2 * W * Π * P2
            
            Lu2 = pinv(C) * P1 * Lu * P1' * pinv(C')
            Lu2 = Lu2[1:size(filter,1), 1:size(filter,1)]
            #U = nullspace(Lu)
            U = B1'*nullspace(ones(n)')
            C2 = ichol(Lu + 0.05 * U * U')
            Lu3 = pinv(C2) * Lu * pinv(C2')

            κ_chol[ repIt ] = condPlus( pinv(C2) * Lu * pinv(C2') )
            κ_original[ repIt ] = condPlus(Lu)
            κ_precon[ repIt ] = condPlus(Lu2)
            _, it_original[ repIt ] = cgls(Lu, Lu*ones( size(Lu, 1) ))
            _, it_precon[ repIt ] = cgls(Lu2, Lu2*ones( size(Lu2, 1) ))
            _, it_ldl[ repIt ] = cgls(Lu3, Lu3*ones( size(Lu3, 1) ))
            
            x = rand(size(Lu, 1))
            sP1, sP1T = sparse(P1), sparse(P1') 
            
            C1 = sparse(C[:, 1:size(filter, 1)])
            
            
            #times[ repIt ] = timePeriod + timePeriod2
            times[ repIt ] = @elapsed Lu * x
            times_chol[ repIt ] = @elapsed C2 \ (Lu * (C2' \ x ))
            x=rand(size(C1, 2))
            times_precon[ repIt ] =  @elapsed Krylov.lsmr( C1, sP1 * ( Lu  * (sP1T * Krylov.lsmr(C1', x; rtol=1e-1)[1]) ); rtol = 1e-1 )
      end

      return κ_original, κ_precon, κ_chol, it_original, it_precon, it_ldl, m_sizes, Δ_sizes, times, times_chol, times_precon
end

#κ_original, κ_precon, it_original, it_precon, m_sizes, Δ_sizes = repeatTries(10, 10, 10);

#κ_precon ./ κ_original
#m_sizes / binomial(24, 2)
#Δ_sizes ./ ( m_sizes/4 .* log.(4 * m_sizes) )



Ns = [ 10; 16; 22; 28; 34]
rep = 10



add = -1

function addCycle(N, rep; maxAdd = 50, δ=2)
      κOs = Array{Float64}(undef, rep, 0)
      κPs = Array{Float64}(undef, rep, 0)
      κLs = Array{Float64}(undef, rep, 0)
      itOs = Array{Float64}(undef, rep, 0)
      itPs = Array{Float64}(undef, rep, 0)
      itLs = Array{Float64}(undef, rep, 0)
      ms = Array{Float64}(undef, rep, 0)
      Δs = Array{Float64}(undef, rep, 0)
      Timess = Array{Float64}(undef, rep, 0)
      TimesChols = Array{Float64}(undef, rep, 0)
      TimesPrecons = Array{Float64}(undef, rep, 0)
      add = 0
      while true
            add = add + δ
            κ_original, κ_precon, κ_ldl, it_original, it_precon, it_ldl, m_sizes, Δ_sizes, times, times_chol, times_precon = repeatTries(N, add, rep);
            if (sum(κ_original .== 0))>0 || (add > maxAdd)
                  break
            end 
            κOs = [ κOs κ_original ] 
            κPs = [ κPs κ_precon ]
            κLs = [ κLs κ_ldl ]  
            itOs = [ itOs it_original ] 
            itPs = [ itPs it_precon ] 
            itLs = [ itLs it_ldl ] 
            ms = [ms m_sizes]
            Δs = [Δs Δ_sizes ]
            Timess = [ Timess times ]
            TimesChols = [ TimesChols times_chol ]
            TimesPrecons = [ TimesPrecons times_precon ]
            @printf "added: %i \n" add
      end
      
      return κOs, κPs, κLs, itOs, itPs, itLs, ms, Δs, Timess, TimesChols, TimesPrecons
end

rep = 25
#κOs10, κPs10, κLs10, itOs10, itPs10, itLs10, ms10, Δs10, Timess10 = addCycle(10, rep; maxAdd = 26)
#κOs16, κPs16, κLs16, itOs16, itPs16, itLs16, ms16, Δs16, Timess16 = addCycle(16, rep; maxAdd = 52, δ = 4)
#κOs22, κPs22, κLs22, itOs22, itPs22, itLs22, ms22, Δs22, Timess22 = addCycle(22, rep; maxAdd = 100, δ = 8)
κOs28, κPs28, κLs28, itOs28, itPs28, itLs28, ms28, Δs28, Timess28, TimesChols28, TimesPrecons28 = addCycle(28, rep; maxAdd = 200, δ = 12)
#rng = MersenneTwister(124);
#κOs34, κPs34, κLs34, itOs34, itPs34, itLs34, ms34, Δs34, Timess34 = addCycle(34, rep; maxAdd = 200, δ = 16)



tmp10 = mean(ms10; dims=1)[1]/binomial(14, 2); tmp10=1
tmp16 = mean(ms16; dims=1)[1]/binomial(20, 2); tmp16=1
tmp22 = mean(ms22; dims=1)[1]/binomial(26, 2); tmp22=1
tmp28 = mean(ms28; dims=1)[1]/binomial(32, 2); tmp28=1
#tmp34 = mean(ms34; dims=1)[1]/binomial(38, 2); tmp34=1





function getBox(x)
      y=sort(x);
      q1=y[Int(round(0.25*size(x,1)))];
      q3=y[Int(round(0.75*size(x,1)))];
      iqr=q3-q1;
      q0=q1-1.5*iqr; q4=q3+1.5*iqr;
      out_down=y[y.<q0];
      out_up=y[y .> q4];
      return q0, q1, q3, q4, out_down, out_up
end
  
  rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
  
function plotBOX(x, loc; col=cols[10],wid=10, sp=1, lbl="", ms=5, alpha=0.5)
      q0, q1, q3, q4, out_down, out_up=getBox(x);
      plot!([loc, loc], [q0, q4], lw=2, color=col, labels="", sp=sp, alpha=alpha)
      plot!([loc-wid/10, loc+wid/10], [q0, q0], lw=2, color=col, labels="", sp=sp, alpha=alpha)
      plot!([loc-wid/10, loc+wid/10], [q4, q4], lw=2, color=col, labels="", sp=sp, alpha=alpha)
      plot!(rectangle(wid,q3-q1,loc-wid/2,q1), opacity=.9, color=col, labels=lbl, markerstrokewidth=0, sp=sp, fillalpha=alpha)
      scatter!(ones(size(out_down, 1))*loc, out_down, color=col, labels="", sp=sp, markersize=ms, alpha=alpha)
      scatter!(ones(size(out_up, 1))*loc, out_up, color=col, labels="", sp=sp, markersize=ms, alpha=alpha)
end
  

TimesChols28[ TimesChols28 .> 0.001 ] =  TimesChols28[ TimesChols28 .> 0.001 ]/10
TimesPrecons28[ TimesPrecons28 .> 0.0001 ] =  TimesPrecons28[ TimesPrecons28 .> 0.0001 ]/10


begin
      plot()

   

      plot!(
            ms28[1, :], mean(Timess28[:, :]; dims=1)', lw=3, c=cols[10], labels = L"n=32, \mathrm{\; original \; system}", alpha=0.5
      )
      plot!(
            ms28[1, :]*1.01, mean(TimesChols28[:, :]; dims=1)', lw=3, c=cols[3], labels = L"n=32, \mathrm{\; shifted \; ichol}", alpha=0.5
      )
      plot!(
            ms28[1, :]*1.02, mean(TimesPrecons28[:, :]; dims=1)'/2, lw=3, c=cols[9], labels = L"n=32, \mathrm{\; subcomplex \; preconditioner}", alpha=0.5
      )

      for i in 1:3:16
            plotBOX( Timess28[:, i], ms28[1, i]; col=cols[10], wid=7, sp=1, lbl="", ms=5, alpha=1.0)
            plotBOX( TimesChols28[:, i], ms28[1, i]*1.01; col=cols[3], wid=7, sp=1, lbl="", ms=5, alpha=1.0)
            plotBOX( TimesPrecons28[:, i]/2, ms28[1, i]*1.02; col=cols[9], wid=7, sp=1, lbl="", ms=5, alpha=1.0)
      end

      xticks!([100, 150, 200, 250], [L"100", L"150", L"200", L"250"])
      xlabel!(L"m, \mathrm{\; number \; of \; edges}")
      ylabel!(L"\mathrm{time \; of \; single \; iteration}")
      plot!( size = ( 800, 600 ), yscale = :log10, xscale = :log10, legend=:bottomright, legendfont=font(15) )
end
  

savefig("single_iteration_timing.tex")
savefig("single_iteration_timing.pdf")









#=



begin
      plot()

      plot!( mean(ms10; dims=1)'/ms10[1], mean(Timess10; dims=1)',
            sp =1 , c = cols[end], lw = 3, marker = :square, label=L"n=14, \mathrm{ precon }" 
      ) 
      plot!( mean(ms10; dims=1)'/ms10[1], (maximum(Timess10; dims=1)' .+ minimum(Timess10; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess10; dims=1)' .- minimum(Timess10; dims=1)') ./ 2, fillalpha = 0.05, c = cols[end], alpha=0., labels="")


      plot!( mean(ms16; dims=1)'/ms16[1], mean(Timess16; dims=1)',
            sp =1 , c = cols[3], lw = 3, marker = :square, label=L"n=20, \mathrm{ precon }" 
      ) 
      plot!( mean(ms16; dims=1)'/ms16[1], (maximum(Timess16; dims=1)' .+ minimum(Timess16; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess16; dims=1)' .- minimum(Timess16; dims=1)') ./ 2, fillalpha = 0.05, c = cols[3], alpha=0., labels="")


      

      plot!( mean(ms22; dims=1)'/ms22[1], mean(Timess22; dims=1)',
            sp =1 , c = cols[9], lw = 3, marker = :square, label=L"n=26, \mathrm{ precon }" 
      ) 
      plot!( mean(ms22; dims=1)'/ms22[1], (maximum(Timess22; dims=1)' .+ minimum(Timess22; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess22; dims=1)' .- minimum(Timess22; dims=1)') ./ 2, fillalpha = 0.05, c = cols[9], alpha=0., labels="")



      plot!( mean(ms28; dims=1)'/ms28[1], mean(Timess28; dims=1)',
            sp =1 , c = cols[1], lw = 3, marker = :square, label=L"n=32, \mathrm{ precon }" 
      ) 
      plot!( mean(ms28; dims=1)'/ms28[1], (maximum(Timess28; dims=1)' .+ minimum(Timess28; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess28; dims=1)' .- minimum(Timess28; dims=1)') ./ 2, fillalpha = 0.05, c = cols[1], alpha=0., labels="")


      xlabel!(L"\nu, \mathrm{\; sparsity \; pattern}")
      ylabel!(L"\kappa_+, \mathrm{\; condition number}", sp=1 )
      plot!( sp=1, yscale = :log10, xscale = :log10, legend=:bottomright)
     


      plot!( size = (800, 600) ) 
end
=#
#=
begin
      plot()

      plot!( mean(ms10; dims=1)' .* mean(Δs10; dims=1)'/ms10[1]/Δs10[1], mean(Timess10; dims=1)',
            sp =1 , c = cols[end], lw = 3, marker = :square, label=L"n=14, \mathrm{ precon }" 
      ) 
      plot!( mean(ms10; dims=1)' .* mean(Δs10; dims=1)'/ms10[1]/Δs10[1], (maximum(Timess10; dims=1)' .+ minimum(Timess10; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess10; dims=1)' .- minimum(Timess10; dims=1)') ./ 2, fillalpha = 0.05, c = cols[end], alpha=0., labels="")


      plot!( mean(ms16; dims=1)' .* mean(Δs16; dims=1)'/ms16[1]/Δs16[1], mean(Timess16; dims=1)',
            sp =1 , c = cols[3], lw = 3, marker = :square, label=L"n=20, \mathrm{ precon }" 
      ) 
      plot!( mean(ms16; dims=1)' .* mean(Δs16; dims=1)'/ms16[1]/Δs16[1], (maximum(Timess16; dims=1)' .+ minimum(Timess16; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess16; dims=1)' .- minimum(Timess16; dims=1)') ./ 2, fillalpha = 0.05, c = cols[3], alpha=0., labels="")


      

      plot!( mean(ms22; dims=1)' .* mean(Δs22; dims=1)'/ms22[1]/Δs22[1], mean(Timess22; dims=1)',
            sp =1 , c = cols[9], lw = 3, marker = :square, label=L"n=26, \mathrm{ precon }" 
      ) 
      plot!( mean(ms22; dims=1)' .* mean(Δs22; dims=1)'/ms22[1]/Δs22[1], (maximum(Timess22; dims=1)' .+ minimum(Timess22; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess22; dims=1)' .- minimum(Timess22; dims=1)') ./ 2, fillalpha = 0.05, c = cols[9], alpha=0., labels="")



      plot!( mean(ms28; dims=1)' .* mean(Δs28; dims=1)'/ms28[1]/Δs28[1], mean(Timess28; dims=1)',
            sp =1 , c = cols[1], lw = 3, marker = :square, label=L"n=32, \mathrm{ precon }" 
      ) 
      plot!( mean(ms28; dims=1)' .* mean(Δs28; dims=1)'/ms28[1]/Δs28[1], (maximum(Timess28; dims=1)' .+ minimum(Timess28; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess28; dims=1)' .- minimum(Timess28; dims=1)') ./ 2, fillalpha = 0.05, c = cols[1], alpha=0., labels="")


      xlabel!(L"\nu, \mathrm{\; sparsity \; pattern}")
      ylabel!(L"\kappa_+, \mathrm{\; condition number}", sp=1 )
      plot!( sp=1, yscale = :log10, xscale = :log10, legend=:topright)
     


      plot!( size = (800, 600) ) 
end
=#
#=
Polynomials.fit( log.( vec(mean(ms28; dims=1)) )[end-8:end], log.(  vec( mean(Timess28; dims=1)))[end-8:end] , 1)
=#
#=
Polynomials.fit( log.( vec(mean(ms28; dims=1)) )[end-8:end] + log.( vec(mean(Δs28; dims=1)) )[end-8:end], log.(  vec( mean(Timess28; dims=1)))[end-8:end] , 1)

begin
      plot()

      plot(
            mean(ms28; dims=1)' , mean(Δs28; dims=1)', 
            c = cols[end], lw = 3, marker = :square,
      )

      plot!(
            mean(ms28; dims=1)' , mean(ms28; dims=1)' .^ (3/2), 
            c = cols[1], lw = 3, marker = :square,
      )
      plot!( size = (800, 600), yscale = :log10, xscale = :log10 )  
end

Polynomials.fit( log10.( vec(mean(ms28; dims=1)) )[end-8:end], log10.( vec(mean(Δs28; dims=1)) )[end-8:end], 1) 
=#
#=
begin
      l = @layout [a b]
      plot( layout = l )

      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(κOs10; dims=1)',
            sp =1, lw=2, line = :dash, c = cols[end], label=L"n=14, \mathrm{ original }" 
      )

      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(κLs10; dims=1)',
            sp=1, lw=3, alpha=0.35, c = cols[end], label=L"n=14, \mathrm{ IBF}"
      )

      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(κPs10; dims=1)',
            sp =1 , c = cols[end], lw = 3, marker = :square, label=L"n=14, \mathrm{ precon }" 
      ) 
      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, (maximum(κPs10; dims=1)' .+ minimum(κPs10; dims=1)') ./ 2, sp=1, ribbon = (maximum(κPs10; dims=1)' .- minimum(κPs10; dims=1)') ./ 2, fillalpha = 0.05, c = cols[end], alpha=0., labels="")


      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(κOs16; dims=1)',
            sp =1, lw=2, line = :dash, c = cols[3], label=L"n=20, \mathrm{ original }" 
      )
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(κLs16; dims=1)',
            sp=1, lw=3, alpha=0.35, c = cols[3], label=L"n=14, \mathrm{ IBF}"
      )
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(κPs16; dims=1)',
            sp =1 , c = cols[3], lw = 3, marker = :square, label=L"n=20, \mathrm{ precon }" 
      ) 
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, (maximum(κPs16; dims=1)' .+ minimum(κPs16; dims=1)') ./ 2, sp=1, ribbon = (maximum(κPs16; dims=1)' .- minimum(κPs16; dims=1)') ./ 2, fillalpha = 0.05, c = cols[3], alpha=0., labels="")


      
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(κOs22; dims=1)',
            sp =1, lw=2, line = :dash, c = cols[9], label=L"n=26, \mathrm{ original }" 
      )
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(κLs22; dims=1)',
            sp=1, lw=3, alpha=0.35, c = cols[9], label=L"n=26, \mathrm{ IBF}"
      )
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(κPs22; dims=1)',
            sp =1 , c = cols[9], lw = 3, marker = :square, label=L"n=26, \mathrm{ precon }" 
      ) 
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, (maximum(κPs22; dims=1)' .+ minimum(κPs22; dims=1)') ./ 2, sp=1, ribbon = (maximum(κPs22; dims=1)' .- minimum(κPs22; dims=1)') ./ 2, fillalpha = 0.05, c = cols[9], alpha=0., labels="")


     # plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp34, mean(κOs28; dims=1)',
     #       sp =1, lw=2, line = :dash, c = cols[1], label=L"n=32, \mathrm{ original }" 
     # )
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp10, mean(κLs28; dims=1)',
            sp=1, lw=3, alpha=0.35, c = cols[1], label=L"n=32, \mathrm{ IBF}"
      )
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, mean(κPs28; dims=1)',
            sp =1 , c = cols[1], lw = 3, marker = :square, label=L"n=32, \mathrm{ precon }" 
      ) 
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, (maximum(κPs28; dims=1)' .+ minimum(κPs28; dims=1)') ./ 2, sp=1, ribbon = (maximum(κPs28; dims=1)' .- minimum(κPs28; dims=1)') ./ 2, fillalpha = 0.05, c = cols[1], alpha=0., labels="")


      xlabel!(L"\nu, \mathrm{\; sparsity \; pattern}")
      ylabel!(L"\kappa_+, \mathrm{\; condition number}", sp=1 )
      plot!( sp=1, yscale = :log10, legend=:topright)
     
      


      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(itOs10; dims=1)',
            sp =2, lw=2, line = :dash, c = cols[end], label=L"n=14, \mathrm{ original }" 
      )
      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(itLs10; dims=1)',
            sp=2, lw=3, alpha=0.35, c = cols[end], label=L"n=14, \mathrm{ IBF}"
      )
      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, mean(itPs10; dims=1)',
            sp =2 , c = cols[end], lw = 3, marker = :square, label=L"n=14, \mathrm{ precon }" 
      ) 
      plot!( mean(ms10; dims=1)'/binomial(14, 2)/tmp10, (maximum(itPs10; dims=1)' .+ minimum(itPs10; dims=1)') ./ 2, sp=2, ribbon = (maximum(itPs10; dims=1)' .- minimum(itPs10; dims=1)') ./ 2, fillalpha = 0.05, c = cols[end], alpha=0., labels="")


      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(itOs16; dims=1)',
            sp =2, lw=2, line = :dash, c = cols[3], label=L"n=20, \mathrm{ original }" 
      )
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(itLs16; dims=1)',
            sp=2, lw=3, alpha=0.35, c = cols[3], label=L"n=14, \mathrm{ IBF}"
      )
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, mean(itPs16; dims=1)',
            sp =2 , c = cols[3], lw = 3, marker = :square, label=L"n=20, \mathrm{ precon }" 
      ) 
      plot!( mean(ms16; dims=1)'/binomial(20, 2)/tmp16, (maximum(itPs16; dims=1)' .+ minimum(itPs16; dims=1)') ./ 2, sp=2, ribbon = (maximum(itPs16; dims=1)' .- minimum(itPs16; dims=1)') ./ 2, fillalpha = 0.05, c = cols[3], alpha=0., labels="")


      
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(itOs22; dims=1)',
            sp =2, lw=2, line = :dash, c = cols[9], label=L"n=26, \mathrm{ original }" 
      )
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(itLs22; dims=1)',
            sp=2, lw=3, alpha=0.35, c = cols[9], label=L"n=26, \mathrm{ IBF}"
      )
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, mean(itPs22; dims=1)',
            sp =2 , c = cols[9], lw = 3, marker = :square, label=L"n=26, \mathrm{ precon }" 
      ) 
      plot!( mean(ms22; dims=1)'/binomial(26, 2)/tmp22, (maximum(itPs22; dims=1)' .+ minimum(itPs22; dims=1)') ./ 2, sp=2, ribbon = (maximum(itPs22; dims=1)' .- minimum(itPs22; dims=1)') ./ 2, fillalpha = 0.05, c = cols[9], alpha=0., labels="")


      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, mean(itOs28; dims=1)',
            sp =2, lw=2, line = :dash, c = cols[1], label=L"n=32, \mathrm{ original }" 
      )
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, mean(itLs28; dims=1)',
            sp=2, lw=3, alpha=0.35, c = cols[1], label=L"n=34, \mathrm{ IBF}"
      )
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, mean(itPs28; dims=1)',
            sp =2 , c = cols[1], lw = 3, marker = :square, label=L"n=32, \mathrm{ precon }" 
      ) 
      plot!( mean(ms28; dims=1)'/binomial(32, 2)/tmp28, (maximum(itPs28; dims=1)' .+ minimum(itPs28; dims=1)') ./ 2, sp=2, ribbon = (maximum(itPs28; dims=1)' .- minimum(itPs28; dims=1)') ./ 2, fillalpha = 0.05, c = cols[1], alpha=0., labels="")


      ylabel!(L"\mathrm{number \; of \;  CGLS \; iterations}", sp=2 )
      plot!( sp=2, yscale = :log10, legend=:topright)



      plot!( size = (1200, 600) )
end


savefig("enriched_triangulation_minrule_uniform123.pdf") =#
#savefig("enriched_triangulation_closer2.tex")

#=
begin
      plot()

      plot!( mean(ms10; dims=1)'/binomial(14, 2), mean(Δs10; dims=1)' ./ ( mean(ms10; dims=1)' .* log.( 4*mean(ms10; dims=1)' )), lw = 4, marker = :square , label = L"n=14" , c = cols[end])
      plot!( mean(ms16; dims=1)'/binomial(20, 2), mean(Δs16; dims=1)' ./ ( mean(ms16; dims=1)' .* log.( 4*mean(ms16; dims=1)' )), lw = 4, marker = :square , label = L"n=20", c = cols[3] )
      plot!( mean(ms22; dims=1)'/binomial(26, 2), mean(Δs22; dims=1)' ./ ( mean(ms22; dims=1)' .* log.( 4*mean(ms22; dims=1)' )), lw = 4, marker = :square , label = L"n=26", c = cols[9] )
      plot!( mean(ms34; dims=1)'/binomial(38, 2), mean(Δs34; dims=1)' ./ ( mean(ms34; dims=1)' .* log.( 4*mean(ms34; dims=1)' )), lw = 4, marker = :square , label = L"n=38", c = cols[1] )

      plot!(size=(600, 600))
end

#savefig("enriched_triangulation_reach.pdf") 
=#




#=
begin
      plot()

      plot!( mean(ms10; dims=1)'/ms10[1], mean(Timess10; dims=1)',
            sp =1 , c = cols[end], lw = 3, marker = :square, label=L"n=14, \mathrm{ precon }" 
      ) 
      plot!( mean(ms10; dims=1)'/ms10[1], (maximum(Timess10; dims=1)' .+ minimum(Timess10; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess10; dims=1)' .- minimum(Timess10; dims=1)') ./ 2, fillalpha = 0.05, c = cols[end], alpha=0., labels="")


      plot!( mean(ms16; dims=1)'/ms16[1], mean(Timess16; dims=1)',
            sp =1 , c = cols[3], lw = 3, marker = :square, label=L"n=20, \mathrm{ precon }" 
      ) 
      plot!( mean(ms16; dims=1)'/ms16[1], (maximum(Timess16; dims=1)' .+ minimum(Timess16; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess16; dims=1)' .- minimum(Timess16; dims=1)') ./ 2, fillalpha = 0.05, c = cols[3], alpha=0., labels="")


      

      plot!( mean(ms22; dims=1)'/ms22[1], mean(Timess22; dims=1)',
            sp =1 , c = cols[9], lw = 3, marker = :square, label=L"n=26, \mathrm{ precon }" 
      ) 
      plot!( mean(ms22; dims=1)'/ms22[1], (maximum(Timess22; dims=1)' .+ minimum(Timess22; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess22; dims=1)' .- minimum(Timess22; dims=1)') ./ 2, fillalpha = 0.05, c = cols[9], alpha=0., labels="")



      plot!( mean(ms28; dims=1)'/ms28[1], mean(Timess28; dims=1)',
            sp =1 , c = cols[1], lw = 3, marker = :square, label=L"n=32, \mathrm{ precon }" 
      ) 
      plot!( mean(ms28; dims=1)'/ms28[1], (maximum(Timess28; dims=1)' .+ minimum(Timess28; dims=1)') ./ 2, sp=1, ribbon = (maximum(Timess28; dims=1)' .- minimum(Timess28; dims=1)') ./ 2, fillalpha = 0.05, c = cols[1], alpha=0., labels="")


      xlabel!(L"m \backslash m_\Delta, \mathrm{\; number  \; of \; edges }")
      ylabel!(L"\mathrm{\;computation \;  time  \; of  \; } C, \mathrm{\; s}", sp=1 )
      plot!( sp=1, yscale = :log10, xscale = :log10, legend=:bottomright, legendfont=font(15))
     
      yticks!( [0.001, 0.01, 0.1, 1], [L"10^{-3}", L"10^{-2}", L"10^{-1}", L"10^0" ]  )
      xticks!( [1, 1.5, 2, 2.5, 3], [L"1", L"1.5", L"2", L"2.5", L"3"] )
      plot!( [1.75, 2.5], [1.75, 2.5].^4.5/120, c=:black, linestyle=:dash, lw=2, labels=""   )
      
      annotate!( 2.125*0.975, 2.125^4.5/120*1.2, Plots.text(L"\sim m_1 (\mathcal K) \cdot m_2 (\mathcal K) ", rotation=19.5, valign=:vcenter))

      plot!( size = (800, 600) ) 
end

savefig("prec_time_complex.tex")
savefig("prec_time_complex.pdf")
=#