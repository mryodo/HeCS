using LinearAlgebra


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



function ichol(A)
      L = copy(A)
      n = size(A,1)
      for k in  1:n
            #if L[ k, k ] < 0 
            #      L[k, k ] *= -1
            #end
            if abs(L[k,k]) > 1e-8
                  L[k,k] = sqrt(L[k,k])
                  for i = k+1:n
                        if abs( L[i,k] ) > 1e-8
                              L[i,k] = L[i,k]/L[k,k]
                        end
                  end
                  for j in k+1:n
                        for i in j+1:n
                              if abs( L[i,j] ) > 1e-8
                                    L[i,j] = L[i,j]-L[i,k]*L[j,k]
                              end
                        end
                  end
            end
      end
      tril!(L)
      return L
end


L = copy(Lu)
n = size(Lu,1)
k =0

k = k + 1
if abs(L[k,k]) > 1e-8
L[k,k] = sqrt(L[k,k])
for i = k+1:n
    if abs( L[i,k] ) > 1e-8
        L[i,k] = L[i,k]/L[k,k]
    end
end
for j in k+1:n
    for i in j:n
        if abs( L[i,j] ) > 1e-8
            L[i,j] = L[i,j]-L[i,k]*L[j,k]
        end
    end
end
end
diag(L)



 