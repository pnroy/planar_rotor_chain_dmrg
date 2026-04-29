#module expectations
  
#using ITensors
#using LinearAlgebra


#export vN_entropy,polarization,correlation
################################################################################
function vN_entropy(_wf,mbond)
	wf=copy(_wf)	###to ensure that this function is not mutating wf
	if length(wf) == 2
			orthogonalize!(wf, 2)
			U,S,V = svd(wf[2], (siteind(wf,2)))
	else
			orthogonalize!(wf, mbond)
			U,S,V = svd(wf[mbond], (linkind(wf, mbond-1), siteind(wf,mbond)))
	end
	SvN = 0.0
	renyi = 0.0
	schmidtvalues=zeros(dim(S, 1))
	for n=1:dim(S, 1)
			p = S[n,n]^2
			schmidtvalues[n]=p
			SvN -= p * log(p)
		renyi += p^2
	end
	renyi=-0.5*log(renyi)

	return SvN,renyi,schmidtvalues
	end
################################################################################
function polarization(wf,Nsites,Nbasis,evod,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")

if evod == "dvr"
	mux = expect(wf,"X")
	muy = expect(wf,"Y")

	dumX = correlation_matrix(wf,"X","X")
	dumY = correlation_matrix(wf,"Y","Y")
else
	mux = expect(wf,"X")
	#muy = expect(wf,"Ycomp")
	muy = expect(wf,"Y")
	
	dumX = correlation_matrix(wf,"X","X")
	#dumY = correlation_matrix(wf,"Ycomp","Ycomp")
	dumY = correlation_matrix(wf,"Y","Y")
end

mux2=0.0
muy2=0.0
for ii=1:Nsites
	mux2+=dumX[ii,ii]
	if evod == "dvr"
		muy2+=dumY[ii,ii]
	else
		muy2-=dumY[ii,ii]
	end
end

return real(sum(mux)),real(sum(muy))
end
################################################################################
function correlation(wf,Nsites,Nbasis,evod,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")

if evod == "dvr"
	dumX = correlation_matrix(wf,"X","X")
	dumY = correlation_matrix(wf,"Y","Y")
else
	dumX = correlation_matrix(wf,"X","X")
	#dumY = correlation_matrix(wf,"Ycomp","Ycomp")
	#pn hack
	dumY = correlation_matrix(wf,"Y","Y")
end

Xcorr=0.0
Ycorr=0.0
for ii=1:Nsites-1
	Xcorr+=dumX[ii,ii+1]
	if evod == "dvr"
		Ycorr+=dumY[ii,ii+1]
	else
		Ycorr-=dumY[ii,ii+1]
	end

end

return real(Xcorr),real(Ycorr)
end
################################################################################
#end

################################################################################
function binder(psi,Nbasis,Xmat,Ymat)

global X = Xmat
global Y = Ymat
global Nspec = Nbasis

include("operators.jl")


#Build O^2 operators#
Nsites = length(psi)
sites = siteinds(psi)

ampoX = AutoMPO()
ampoY = AutoMPO()

for i=1:Nsites
	fac = (-1.0)^i
	ampoX += fac,"X",i
	ampoY += fac,"Y",i
end

Mx = MPO(ampoX,sites)
My = MPO(ampoY,sites)

#M2*Psi#
MxPsi = replaceprime(contract(Mx, psi), 2 => 1)#apply(Mx,psi)
Mx2Psi = replaceprime(contract(Mx, MxPsi), 2 => 1)#apply(Mx,MxPsi)
MyPsi = replaceprime(contract(My, psi), 2 => 1)#apply(My,psi)
My2Psi = replaceprime(contract(My, MyPsi), 2 => 1)#apply(My,MyPsi)

#Expectation values#
Mx2 = real(inner(psi,Mx2Psi))
Mx4 = real(inner(Mx2Psi,Mx2Psi))

My2 = real(inner(psi,My2Psi))
My4 = real(inner(My2Psi,My2Psi))

binderX = 1.0-Mx4/(3.0*Mx2^2)
binderY = 1.0-My4/(3.0*My2^2)

return Mx2,My2,binderX,binderY
end
################################################################################
