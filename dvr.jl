#module dvr

#export exp_dvr
#############################################################################
function exp_dvr(N)

#Write DVR grid and position operators on grid#
phi = zeros((2*N+1))
X = zeros((2*N+1,2*N+1))
Y = zeros((2*N+1,2*N+1))
V6 = zeros((2*N+1,2*N+1))

for ii=1:(2*N+1)
	phi[ii] = ii*2.0*pi/(2*N+1)
	X[ii,ii] = cos(phi[ii])	
	V6[ii,ii] = cos(6.0*phi[ii])	
	Y[ii,ii] = sin(phi[ii])	
end

#Write kinetic energy matrix in DVR basis#
T = zeros((2*N+1,2*N+1))

for ii=1:(2*N+1)
	for jj=ii:(2*N+1)
		T[ii,jj] = cos(pi*(ii-jj)/(2*N+1))/(2.0*sin(pi*(ii-jj)/(2*N+1))^2)*(-1.)^(ii-jj)
		T[jj,ii] = T[ii,jj]
	end
	T[ii,ii] = N*(N+1)/3
end

return T,X,Y,V6
end
#############################################################################
#end
