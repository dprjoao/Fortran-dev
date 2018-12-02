Program mdir
!********************************************************************************************************************!

!Programa que implementa a RTM com condicao de imagem por correlacao cruzada

!Autor: Joao Pedro dos Prazeres Reis

!Ultima atualizacao: 03/04/2018

!*******************************************************************************************************************!

Implicit None

!------------------------------------------------------Declaraao de variaveis------------------------------------------------------------

Integer::conta_sismo,k,j,i,w,d,c,n,v,Nt,Nx,Nz,conta_snap,comp_byte,nsnap,l,fz,pz,fx,px,prof_reg,Tsismo,Ntsismo
Integer::ntiros,nsismo,f,Nce,Ncd,Esptiro,Nr,offmin1,offmin2,dh,cmin,cmax,pos,contador,ncor,ntcor
Real,Allocatable,dimension(:)::Be,Bd
Real,Parameter :: pi=3.14179265538979d0
Real ,Allocatable, dimension(:,:,:) ::  P, Pfonte, soma
Real ,Allocatable, dimension(:,:) :: mod_vel,sismo,U2,U4,P1,P2,P3,Img,sismo2, Img_pond_src
Real Sx,Sz,source,delta_t,P_zero,t,snap,h,rho,t0,fc,A,f_corte,fator_cerjan,fator_c,Lcab
Real :: c1,c2, dummy
character(len=3) :: csnap, ans
character(len=25):: modelo_vp,numsismo,modelo_vps

!-----------------------------------------Inteiro referente ao numero de bytes de um real simples------------------------------------------

comp_byte=1

!-----------------------------------------Começa a contar o tempo de processamento do programa----------------------------------------------

call cpu_time(c1)

!--------------------------------------------------Leitura do arquivo de entrada------------------------------------------------------------

open(12,file='Entrada.txt')
read(12,*) Nx,Nz !Tamanho do grid em x,z 
read(12,*) dh !Espaçamento entre tiros
read(12,*) Ntiros, Lcab, Nr !Numero de tiros, tamanho do cabo e numero de canais
read(12,*) t,Tsismo         !tempo de modelagem
read(12,*) snap
read(12,*) f_corte,A,fx,fz,prof_reg !frequencia de corte, amplitude e posiçao da fonte e do poco
read(12,*) Nce,Ncd,fator_c !Posiçao camada de cerjan em metros e fator de cerjan
read(12,*) modelo_vps,modelo_vp !modelo de velocidades
read(12,*) ncor !numero de amostras do campo de ondas

!----------------------------------Alocacao dinamica do array para o modelo de velocidades--------------------------------------

Allocate(mod_vel(Nz,Nx))

!---------------------------------------Leitura do arquivo de modelo de velocidades--------------------------------

open(100,file=modelo_vps,status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
read(100,rec=1) ((mod_vel(i,j),i=1,Nz),j=1,Nx)
close(100)


!------------------------------------------------Paramentros da modelagem----------------------------------------------------

cmin = MINVAL(mod_vel)
cmax = MAXVAL(mod_vel)
print*, cmax, cmin, modelo_vps
fator_cerjan = fator_c/10000
h = nint(cmin/(f_corte*4)-1)
!dh = nint(esptiro/h)
delta_t = h/(cmax*5)
nt = nint(t/delta_t)
ntcor = nint((t/ncor)/delta_t)
nsnap = nint((t/snap)/delta_t)
t0 = 2*sqrt(pi)/f_corte
fc = f_corte/(3*sqrt(pi))
Px = Nx*h
Pz = Nz*h
offmin2 = Lcab/Nr
Ntsismo = nt/Tsismo
dummy = 0.00000001

!----------------------------------Alocacao dinamica de variaveis--------------------------------------

Allocate (P1(Nz,Nx),P2(Nz,Nx),P3(Nz,Nx),U2(Nz,Nx),U4(Nz,Nx),Be(Nce),Bd(Ncd) &
& ,P(Nz,Nx,Ntcor),Img(Nz,Nx),sismo(Ntsismo,Nr),sismo2(Ntsismo,Nz),Pfonte(Nz,Nx,Ntcor),soma(Nz,Nx,Ntcor),Img_pond_src(Nz,Nx))

!----------------------------------Vizualizaçao dos parametros da modelagem--------------------------------------
print*,'Tamanho do modelo(m) e tempo total(s) ---->', Px, Pz, t
pause

print*,'Espacamento entre pontos e dt ---->', h, delta_t
pause

print*,'Passos de tempo, Passos de tempo array de correlacao cruzada ---->', Nt, Ntcor
pause

print*,'Numero de  tiros, offsetminimo (m), tamanho de cabo (m) e numero de canais ---->', ntiros, offmin2*h, Lcab, Nr
pause

print*,'Continuar com a RTM?[s/n]'
read(*,*) ans

if ((ans == 'n')) then
   print*, 'RTM cancelada!'
   stop   
endif
!---------------------------------Vetorizacao do amortecimentos(Cerjan)------------------------------------

Be=0
do c =1,Nce
Be(c)=exp(-(fator_cerjan*(Nce-c))*(fator_cerjan*(Nce-c)))
enddo
Bd=0
do c =1,Ncd
Bd(c)=exp(-(fator_cerjan*(Ncd-c))*(fator_cerjan*(Ncd-c)))
enddo


!-----------------------------------------Segunda ordem------------------------------------------------

U2=0
do j=2,Nx-1
   do i=2,Nz-1
	  U2(i,j)=(mod_vel(i,j)**2)*delta_t**2/(h*h)
   enddo
enddo

!------------------------------------------Quarta ordem------------------------------------------------

U4=0
do j=3,Nx-2
   do i=3,Nz-2
   U4(i,j)=(mod_vel(i,j)**2)*delta_t**2/(12*h*h)
   enddo
enddo

!------------------------------------------Inicio do loop dos tiros------------------------------------

conta_sismo=0
contador=0.0
img=0.0
nsismo=0
img_pond_src = 0.0
soma = 0.0

do f=1,ntiros
   conta_sismo=conta_sismo+1
   conta_snap=0

!----------------------------------Chamada de rotina para as condicoes iniciais------------------------------------------   
	sismo =0.0
	P = 0.0
	Pfonte = 0.0
   call icond()  

!-------------------------------------------Abrindo sismogramas---------------------------------------

   write(numsismo,"(I3.3)")f
   Print*, "Abrindo sismo",f
   open(11,file='sismo'//trim(adjustl(numsismo))//'.bin',status="old",access="direct",form="unformatted",recl=(comp_byte*Nt*Nr))
   read(11,rec=1) ((sismo(i,j),i=1,Nt),j=1,Nr)
   close(11)
   !open(7,File="vsp"//trim(adjustl(numsismo))//".bin",status="replace",access="direct",form="unformatted",recl=(comp_byte*Ntsismo*Nz))
   !write(7,rec=1) ((sismo2(i,j),i=1,Ntsismo),j=1,Nz)

!-----------------------------------------Inicico do loop do tempo------------------------------------

   do k=1,Nt
   !	if (mod(k,100)==0) then
   !       write(*,*) k*delta_t ,  'segundos'
    !  endif

!-------------------------------------------Calculo do termo fonte----------------------------------------------------

      source= A*(2*pi*(pi*fc*(k*delta_t-t0))*(pi*fc*(k*delta_t-t0))-1.0)*exp(-pi*(pi*fc*(k*delta_t-t0))**2)

!-------------------------------------------Injecao do pulso sismico no modelo-----------------------------------------------------      

      P1(fz,fx+(dh*(f-1)))=P1(fz,fx+(dh*(f-1)))+source
     
!--------------------------------------------Calculo dos paineis de pressao--------------------------------------------------

      !diferenças finitas de segunda ordem nas bordas

      do j=2,Nx-1 !borda superior
	 		P3(2,j)= U2(2,j)*(P2(3,j)+P2(1,j)-4*P2(2,j)+P2(2,j+1)+P2(2,j-1))+2*P2(2,j)-P1(2,j)
      enddo

      do j=2,Nx-1 !borda inferior
	 		P3(nz-1,j)= U2(nz-1,j)*(P2(nz,j)+P2(nz-2,j)-4*P2(nz-1,j)+P2(nz-1,j+1)+P2(nz-1,j-1))+2*P2(nz-1,j)-P1(nz-1,j)
      enddo

      do i=2,Nz-1 !borda a esquerda
	 		P3(i,2)= U2(i,2)*(P2(i+1,2)+P2(i-1,2)-4*P2(i,2)+P2(i,3)+P2(i,1))+2*P2(i,2)-P1(i,2)
      enddo

      do i=2,Nz-1 !borda a direita
	 		P3(i,nx-1)= U2(i,nx-1)*(P2(i+1,nx-1)+P2(i-1,nx-1)-4*P2(i,nx-1)+P2(i,nx)+P2(i,nx-2))+2*P2(i,nx-1)-P1(i,nx-1)
      enddo

      do j=3,Nx-2
         do i=3,Nz-2
            P3(i,j)= U4(i,j)*(-P2(i+2,j)-P2(i-2,j)-P2(i,j+2)-P2(i,j-2)-60*P2(i,j)+16*(P2(i-1,j)+P2(i,j-1)+P2(i+1,j)&
            +P2(i,j+1)))+2*P2(i,j)-P1(i,j)
         enddo
      enddo

      !campo de pressoes acusticas
      do w=1,Nce-1 !esquerda
         do d=1,Nz
            P3(d,w)=Be(w)*P3(d,w)
            P2(d,w)=Be(w)*P2(d,w)
         enddo
      enddo

      n=Ncd
      do w=Nx-Ncd+2,Nx  !direita
         n=n-1
         do d=1,Nz
            P3(d,w)=Bd(n)*P3(d,w)
            P2(d,w)=Bd(n)*P2(d,w)
         enddo
      enddo
      
      do d=1,Nx
         n=Nce
         do w=Nz-Nce+2,Nz !inferior
            n=n-1
            P3(w,d)=Be(n)*P3(w,d)
            P2(w,d)=Be(n)*P2(w,d)
         enddo
      enddo

      !Reynolds
      do n=1,nz !Nao reflexao da borda esquerda
   	     P3(n,1)=(delta_t*mod_vel(n,1)/h)*(P2(n,2)-P2(n,1))+P2(n,1)
      enddo

      do n=1,nz !Nao reflexao da borda direita
   	     P3(n,nx)= -(delta_t*mod_vel(n,nx)/h)*(P2(n,nx)-P2(n,nx-1))+P2(n,nx)
      enddo

      do i=2,nx-1 !Nao reflexao da borda debaixo
   	     P3(nz,i)= -(delta_t*mod_vel(nz,i)/h)*(P2(nz,i)-P2(nz-1,i))+P2(nz,i)
      enddo

      P3(1,:)=0.0!Superficie livre
     
    !  if (mod(k,nsnap)==0) then
   !	   conta_snap=conta_snap+1
     !    write(csnap,"(I2.1)")conta_snap
        !Abrindo o arquivo de saida
     !    open(6,File=trim(csnap)//".bin",status="replace",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
         !imprimindo no arquivo de saida os resultados
     !    write(6,rec=1) ((P3(i,j),i=1,Nz),j=1,Nx)
     !endif
      
      !Costruindo a matriz de correlacao cruzada
      if (mod(k,ncor)==0) then
          contador=contador+1
          do j=1,Nx
      	    do i=1,Nz
                P(i,j,contador)=P3(i,j) 
					 Pfonte(i,j,contador) = P(i,j,contador)**2
             enddo
          enddo
      endif
      
      !Marcha no tempo e atualizacao dos campos
      P1=P2
      P2=P3
      

!*------------------------------------------------------fim do loop de tempo-----------------------------------------------------------*!

	enddo
   !write(*,*) "teste"
   !open(22,file='teste.bin',status="unknown",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
   !Do k=1,10  
!	write(22,rec=k) ((P(i,j,100*k),i=1,Nz),j=1,Nx)
!   enddo

!----------------------------------------------Chamada de rotina para as condicoes iniciais------------------------------------------   

   call icond()

!-------------------------------------------------------Inicio do loop do tempo----------------------------------------------------------
	 do k=Nt,1,-1
		! if (mod(k,60)==0) then
		 !   print*, k*delta_t,  'segundos'
		 !endif
		 
       !do i=2,Nz
		    
        !  P1(i,pos)=P1(i,pos)+sismo2(k,i)
		 
       !enddo   
		 
       do j=1,76

		    P1(prof_reg,(fx-(j*offmin2))+((f-1)*dh))=P1(prof_reg,(fx-(j*offmin2))+((f-1)*dh))+sismo(k,j)

		 enddo
!--------------------------------------------Calculo dos paineis de pressao--------------------------------------------------

      !diferenças finitas de segunda ordem nas bordas

      do j=2,Nx-1 !borda superior
	 		P3(2,j)= U2(2,j)*(P2(3,j)+P2(1,j)-4*P2(2,j)+P2(2,j+1)+P2(2,j-1))+2*P2(2,j)-P1(2,j)
      enddo

      do j=2,Nx-1 !borda inferior
	 		P3(nz-1,j)= U2(nz-1,j)*(P2(nz,j)+P2(nz-2,j)-4*P2(nz-1,j)+P2(nz-1,j+1)+P2(nz-1,j-1))+2*P2(nz-1,j)-P1(nz-1,j)
      enddo

      do i=2,Nz-1 !borda a esquerda
	 		P3(i,2)= U2(i,2)*(P2(i+1,2)+P2(i-1,2)-4*P2(i,2)+P2(i,3)+P2(i,1))+2*P2(i,2)-P1(i,2)
      enddo

      do i=2,Nz-1 !borda a direita
	 		P3(i,nx-1)= U2(i,nx-1)*(P2(i+1,nx-1)+P2(i-1,nx-1)-4*P2(i,nx-1)+P2(i,nx)+P2(i,nx-2))+2*P2(i,nx-1)-P1(i,nx-1)
      enddo

      do j=3,Nx-2
         do i=3,Nz-2
            P3(i,j)= U4(i,j)*(-P2(i+2,j)-P2(i-2,j)-P2(i,j+2)-P2(i,j-2)-60*P2(i,j)+16*(P2(i-1,j)+P2(i,j-1)+P2(i+1,j)&
            +P2(i,j+1)))+2*P2(i,j)-P1(i,j)
         enddo
      enddo

      !campo de pressoes acusticas
      do w=1,Nce-1 !esquerda
         do d=1,Nz
            P3(d,w)=Be(w)*P3(d,w)
            P2(d,w)=Be(w)*P2(d,w)
         enddo
      enddo

      n=Ncd
      do w=Nx-Ncd+2,Nx  !direita
         n=n-1
         do d=1,Nz
            P3(d,w)=Bd(n)*P3(d,w)
            P2(d,w)=Bd(n)*P2(d,w)
         enddo
      enddo
      
      do d=1,Nx
         n=Nce
         do w=Nz-Nce+2,Nz !inferior
            n=n-1
            P3(w,d)=Be(n)*P3(w,d)
            P2(w,d)=Be(n)*P2(w,d)
         enddo
      enddo

      !Reynolds
      do n=1,nz !Nao reflexao da borda esquerda
   	     P3(n,1)=(delta_t*mod_vel(n,1)/h)*(P2(n,2)-P2(n,1))+P2(n,1)
      enddo

      do n=1,nz !Nao reflexao da borda direita
   	     P3(n,nx)= -(delta_t*mod_vel(n,nx)/h)*(P2(n,nx)-P2(n,nx-1))+P2(n,nx)
      enddo

      do i=2,nx-1 !Nao reflexao da borda debaixo
   	     P3(nz,i)= -(delta_t*mod_vel(nz,i)/h)*(P2(nz,i)-P2(nz-1,i))+P2(nz,i)
      enddo

      P3(1,:)=0.0!Superficie livre	

	  !Abrindo o arquivo de saida

		  !if (mod(k,nsnap)==0) then
		   !  conta_snap=conta_snap+1
		    ! write(csnap,"(I3.3)")conta_snap

		    !Abrindo o arquivo de saida
		 
           !open(6,File="d"//trim(adjustl(csnap))//".bin",status="replace",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))
		 
          !imprimindo no arquivo de saida os resultados
		 
           !write(6,rec=1) ((P3(i,j),i=1,Nz),j=1,Nx)
		     !close(6)    
		  !endif 

		  if (mod(k,ncor)==0) then
		     do j=1,Nx
		        do i=1,Nz
		           Img(i,j)=(Img(i,j)+P(i,j,contador)*P3(i,j))
					  !soma(i,j,contador) = soma(i,j,contador) + Pfonte(i,j,contador)
                 !Img_pond_src(i,j) = Img(i,j)/(soma(i,j,contador) + dummy)
		        enddo
             enddo
				 contador=contador-1 
		  endif

		  !sobrescrever a matriz
		  P1=P2
		  P2=P3
	  enddo

!-----------------------------------------------imprimindo imagem---------------------------------------------------------!
     print*, "Acumulando imagem", f
     open(8,File="imagem.bin",status="replace",access="direct",form="unformatted",recl=(comp_byte*Nx*Nz))!//trim(adjustl(numsismo))//"
	  write(8,rec=1) ((Img(i,j),i=1,Nz),j=1,Nx)
	  close(8)

!------------------------------------------------Fim do loop de tiros----------------------------------------------------!

enddo

!--------------------------------------------------------------------------------------------------------------------!

!dealocacao de variaveis
deallocate(P1,P2,P3,mod_vel,U2,U4,Be,Bd,P,Pfonte,soma,Img,Img_pond_src)

call cpu_time(c2)
print*, "Tempo total de processamento = ",c2-c1
pause

!----------------------------------------------Subrotina para condicoes iniciais-----------------------------------------
Contains
SUBROUTINE icond()
Implicit None

P1=0.0
P2=0.0
P3=0.0
source=0.0


end Subroutine icond

!-----------------------------------------------------Fim do programa--------------------------------------------------------

end program mdir
