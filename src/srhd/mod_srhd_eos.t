!Adapted for SRHD, D. Millas, January 2019
!Original for SRMHD, Z. Meliani, simplified by R. Keppens
module mod_srhd_eos
  use mod_global_parameters
  use mod_srhd_parameters
  implicit none

 contains

!!DM Modify without the BField 

  !>compute the enthalpy
  subroutine srhd_get_enthalpy(ixO^L,rho,p,rhoh)
    use mod_global_parameters, only: nw, ndim
    integer, intent(in)                :: ixO^L
    double precision, intent(in)       :: rho(ixO^S),p(ixO^S)
    double precision, intent(inout)    :: rhoh(ixO^S)
    double precision, dimension(ixO^S) :: E_th,E

    if(srhd_eos) then
     E_th = p*inv_gamma_1
     E    = E_th+&
                   dsqrt(E_Th**2.0d0+rho**2.0d0)
     rhoh = 0.5d0*((srhd_gamma+1.0d0) * E-&
               gamma_1* rho*(rho/E))
    else
     rhoh = (rho+gamma_to_gamma_1*p)
    end if
  end subroutine srhd_get_enthalpy

  !> Calculate thermal pressure for enthalpy and density
  subroutine srhd_get_pressure_primitive_eos(ixI^L,ixO^L,rho,rhoh,pth)
    use mod_global_parameters
    implicit none
    integer, intent(in)            :: ixI^L, ixO^L
    double precision, intent(in)   :: rho(ixO^S),rhoh(ixO^S)
    double precision, intent(out)  :: pth(ixI^S)

    double precision               :: E(ixO^S)

    cond_iseos : if(srhd_eos) then
     E = (rhoh+dsqrt(rhoh**2.0d0+(srhd_gamma**2.0d0-1.0d0)&
            *rho**2.0d0))/(srhd_gamma+1.0d0)
     pth(ixO^S) = 0.5d0*gamma_1* (E-rho*(rho/E))
    else cond_iseos
     pth(ixO^S) = (rhoh-rho)*inv_gamma_1
    end if cond_iseos
  end subroutine srhd_get_pressure_primitive_eos

  !> Calculate thermal pressure within ixO^L
  subroutine srhd_get_pthermal_eos(ixI^L,ixO^L,x,rho,rhoh,e_in,pth)
    use mod_global_parameters
    implicit none

    integer, intent(in)          :: ixI^L, ixO^L
    double precision, intent(in) :: rho(ixO^S),rhoh(ixO^S),e_in(ixO^S)
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(out):: pth(ixI^S)

    double precision             :: e(ixO^S)

    is_notdiab: if(srhd_energy) then

      is_einternal : if(block%e_is_internal) then
        is_eos_einter : if(srhd_eos)then
         pth(ixO^S) = 0.5d0*gamma_1*(e_in-rho*(rho/e_in))
        else is_eos_einter
         pth(ixO^S)=gamma_1*e_in
        end if is_eos_einter
      else is_einternal
       is_eos : if(srhd_eos) then
         e  = (rhoh+dsqrt(rhoh**2.0d0+(srhd_gamma**2.0d0-1.0d0)&
            *rho**2.0d0))/(srhd_gamma+1.0d0)

         pth(ixO^S) = 0.5d0*gamma_1&
             * (e-rho*(rho/e))
       else is_eos
         pth(ixO^S) = (rhoh - rho)/gamma_to_gamma_1
       end if is_eos

      end if is_einternal

    else is_notdiab

      pth(ixO^S)=srhd_adiab*rho**srhd_gamma

    end if is_notdiab
  end subroutine srhd_get_pthermal_eos


  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  subroutine srhd_get_csound2_eos(ixI^L,ixO^L,x,rho,rhoh,csound2)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: rho(ixO^S),rhoh(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: h(ixO^S),E(ixO^S),p(ixO^S)
    if(srhd_energy) then
      if(srhd_eos) then

       E(ixO^S)=(rhoh(ixO^S)+dsqrt(rhoh(ixO^S)**2.0d0+(srhd_gamma**2.0d0-1.0d0)&
             *rho**2.0d0))/(srhd_gamma+1.0d0)
       p=0.5d0*(E-rho*(rho/E))
       csound2(ixO^S)=(p&
                   *((srhd_gamma+1.0d0)&
                   +gamma_1*(rho/E)**2.0d0))&
                   /(2.0*rhoh)
      else
       csound2(ixO^S)=srhd_gamma*p/rho
      end if
    else
      csound2(ixO^S)=srhd_gamma*srhd_adiab*rho**gamma_1
    end if
  end subroutine srhd_get_csound2_eos

  !> Calculate the square of the thermal sound speed csound2 within ixO^L.
  subroutine srhd_get_csound2_prim_eos(ixI^L,ixO^L,x,rho,rhoh,p,csound2)
    use mod_global_parameters
    implicit none
    integer, intent(in)             :: ixI^L, ixO^L
    double precision, intent(in)    :: rho(ixO^S),rhoh(ixO^S),p(ixO^S)
    double precision, intent(in)    :: x(ixI^S,1:ndim)
    double precision, intent(out)   :: csound2(ixO^S)

    double precision                :: h(ixI^S),E(ixI^S)
    if(srhd_energy) then
      if(srhd_eos) then
       E(ixO^S)=(rhoh(ixO^S)+dsqrt(rhoh(ixO^S)**2.0d0+(srhd_gamma**2.0d0-1.0d0)&
             *rho**2.0d0))/(srhd_gamma+1.0d0)
       csound2(ixO^S)=(p&
                   *((srhd_gamma+1.0d0)&
                   +gamma_1*(rho/E)**2.0d0))&
                   /(2.0*rhoh)
      else
       csound2(ixO^S)=srhd_gamma*p/rho
      end if
    else
      csound2(ixO^S)=srhd_gamma*srhd_adiab*rho**gamma_1
    end if
  end subroutine srhd_get_csound2_prim_eos

  !> Calculate the Enthalpy and dhdp from pressure
  subroutine srhd_get_val_h_dhdp(rho,p,drhodp,h,dhdp)
    use mod_global_parameters
    implicit none

    double precision,           intent(in) :: rho,p,drhodp
    double precision,           intent(out):: h
    double precision, optional, intent(out):: dhdp

    double precision             :: Eth,E,dEthdp,dEdp,sqrtEt2rho2,rhotoE
    
    is_energy : if(srhd_energy) then
      is_eos : if(srhd_eos) then
       
       Eth = p*inv_gamma_1
       sqrtEt2rho2 = dsqrt(Eth**2.0d0+rho**2.0d0) 
       E = (Eth + sqrtEt2rho2)
       rhotoE= rho/E
       h = 0.5 *((srhd_gamma+1.0d0) * E-gamma_1 * rho*rhotoE)

       if(present(dhdp))then
        dEthdp = inv_gamma_1
        dEdp   = dEthdp +(dEthdp*Eth+rho*drhodp)/sqrtEt2rho2
        dhdp   = 0.5 *(((srhd_gamma+1.0d0) + rhotoE**2.0)* dEdp&
                       - 2.0*gamma_1*rhotoE* drhodp)
       end if
      else is_eos
       h=rho+gamma_to_gamma_1*p
       if(present(dhdp))dhdp=drhodp+gamma_to_gamma_1
      end if is_eos
    else is_energy
    end if is_energy
  end  subroutine srhd_get_val_h_dhdp

   !> Calculate the pressure and dpdxi from xi
  subroutine srhd_get_val_p_dpdxi(rho,h,drhodxi,dhdxi,p,dpdxi)
    use mod_global_parameters
    implicit none

    double precision,           intent(in) :: rho,h,drhodxi,dhdxi
    double precision,           intent(out):: p
    double precision,           intent(out):: dpdxi

    double precision                       :: Eth,E,dEdxi,rhotoE

    if(srhd_energy) then
      is_eos : if(srhd_eos) then
       E =(h+dsqrt(h**2.0d0+(srhd_gamma**2.0d0-1.0d0)*rho**2.0d0))&
              /(srhd_gamma+1.0d0)

       rhotoE = rho/E

       ! output pressure
       p=gamma_1/2.0d0* (E-rho*rhotoE)
      
       dEdxi=(dhdxi+(h*dhdxi+(srhd_gamma**2.0d0-1.0d0)*rho*drhodxi)&
             /dsqrt(h**2.0d0+(srhd_gamma**2.0d0-1.0d0)*rho**2.0d0))&
             / (srhd_gamma+1.0d0) 
       
       dpdxi=gamma_1/2.0d0*((1.0+rhotoE**2.0)*dEdxi-2.0*rhotoE*drhodxi)

      else is_eos
       p = (h - rho)/gamma_to_gamma_1
       dpdxi = (dhdxi-drhodxi)/gamma_to_gamma_1

      end if is_eos
    end if 
  end  subroutine srhd_get_val_p_dpdxi

  subroutine srhd_get_h_noflux(rho,h_p,h)
    use mod_global_parameters
    implicit none

    double precision,           intent(in) :: rho,h_p

    double precision,           intent(out):: h

    if(srhd_energy) then
      is_eos : if(srhd_eos) then
       h=0.5*((srhd_gamma+1.0)*h_p-gamma_1*rho*(rho/h_p))  
      else is_eos
       h=rho+(h_p-rho)*srhd_gamma
      end if is_eos
    end if

  end subroutine srhd_get_h_noflux

end module mod_srhd_eos
