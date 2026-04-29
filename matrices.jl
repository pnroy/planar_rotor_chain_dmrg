#module matrices

using LinearAlgebra

#export kinetic,Xoperator,Yoperator,Upoperator,Downoperator,PNoperator,m
#################################################################
function kinetic(mmax)

matrix = zeros((2*mmax+1),(2*mmax+1))

k=0
for m=-mmax:mmax
	k+=1
	matrix[k,k] = m*m
end	

return matrix
end
function m(mmax)
matrix = zeros((2*mmax+1),(2*mmax+1))
k=0
for m=-mmax:mmax
	k+=1
	matrix[k,k] = m
end	
return matrix
end
#################################################################
function Xoperator(mmax)

diagonal = zeros((2*mmax+1))
subdiagonal = 0.5*ones((2*mmax))

matrix=Tridiagonal(subdiagonal,diagonal,subdiagonal)

return matrix
end
#################################################################
function Yoperator(mmax)

diagonal = zeros((2*mmax+1))
upper_subdiagonal = 0.5*ones((2*mmax))
lower_subdiagonal = -0.5*ones((2*mmax))

matrix=Tridiagonal(lower_subdiagonal,diagonal,upper_subdiagonal)

return matrix
end
#################################################################
function Upoperator(mmax)

diagonal = zeros((2*mmax+1))
upper_subdiagonal = zeros((2*mmax))
subdiagonal = ones((2*mmax))
	
matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
	
return matrix
end
#################################################################
function Downoperator(mmax)

diagonal = zeros((2*mmax+1))
upper_subdiagonal = ones((2*mmax))
subdiagonal = zeros((2*mmax))

matrix=Tridiagonal(subdiagonal,diagonal,upper_subdiagonal)
	
return matrix
end

#################################################################


#################################################################
function lz_operator(mmax)

matrix = zeros((2*mmax+1),(2*mmax+1))

k=0
for m=-mmax:mmax
	k+=1
	matrix[k,k] = m
end	

return matrix
end

#################################################################
 


#end 
