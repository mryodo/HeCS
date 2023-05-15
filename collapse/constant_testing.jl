using LinearAlgebra



d = 25
q = 40
M = 0.7*sqrt(d)
EyyT = zeros( d, d )

for i in 1 : d^6
      tmp = 0.7*rand(d, 1)
      EyyT = EyyT + tmp * tmp'
end
EyyT = 1/d^6 * EyyT


rep = q^4

res = 0 

for j in 1:rep
      global d, q, EyyT, res
      Y = 0.7*rand(d, q)


      res = res + opnorm(
            mean( 
                  [
                        Y[:, i] * Y[:, i]' for i in 1:q
                  ]
            ) - EyyT   , 2
      )
end
res = res / rep

res * sqrt(q / log(q )) / M



