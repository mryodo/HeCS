using LinearAlgebra, Arpack


n = 100

x,y = rand(n, 1), rand(n, 1)
x = x/norm(x)
y = y/norm(y)

αs = Array(-1.0:0.01:1.0)
βs = Array(-1.0:0.01:1.0)
Λ = zeros(201, 201);


for i in 1:size(αs, 1)
      global x, y, Λ, αs, βs
      α = αs[i];
      for j in 1:size(βs, 1)
            β = βs[j]
            Mat = (x-y) * (x-y)' + α*x*x' + β*y*y'
            Λ[i, j] = eigs(Mat)[1][1] 
      end
end


using Plots, ColorSchemes, Plots.PlotMeasures,  LaTeXStrings
pgfplotsx()
theme(:mute)
Plots.scalefontsizes(1.25)
cols=ColorSchemes.Spectral_11;

heatmap(Λ)



plot()
plot!( αs, diag(Λ)  )



n = 100

x,y = rand(n, 1), rand(n, 1)
x = x/norm(x)
y = y/norm(y)



αs = Array(1.0:-0.001:-1.0)
λ = zeros(size(αs))


for i in 1:size(αs, 1)
      global x, y, αs, λ
      α = αs[i];
      Mat = -(x-y) * (x-y)' + α*x*x' + α*y*y'
      λ[i] = maximum([abs(eigs(Mat, which=:SR)[1][1]), abs(eigs(Mat, which=:LR)[1][1])])
end

plot()
plot!( αs, abs.(λ), lw=2, )

function δ(i; n = 9)
      res = zeros(n, 1);
      res[i] = 1;
      return res
end


t1 = [ 1; -1; 0; 1; 0; 0; 0; 0; 0 ];
t2 = [ 1; 0; -1; 0; 1; 0; 0; 0; 0 ];

r1 = [ 0; 1; 0; 0; 0; -1; 0; 1; 0 ];
r2 = [ 0; 0; 1; 0; 0; -1; 0; 0; 1 ];
r3 = [ 0; 0; 0; 1; 0; 0; -1; 1; 0 ];
r4 = [ 0; 0; 0; 0; 1; 0; -1; 0; 1 ];

Mat = ( t1 - t2 ) * ( t1 - t2 )'

tmp = - r1 * r1' + r2 * r2' + r3 * r3' - r4 * r4'

( t1 - t2 ) * δ(1) * δ(1)'





