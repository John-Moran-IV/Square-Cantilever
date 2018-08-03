
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi,
     &     sti,xstateini,xstate,nstate_,dtime)
!
!     user subroutine dflux
!
!
!     INPUT:
!
!     sol                current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        1  = body flux
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     temp               currently not used
!     press              currently not used
!     loadtype           load type label
!     area               for surface flux: area covered by the
!                            integration point
!                        for body flux: volume covered by the
!                            integration point
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     co(3,1..nk)        coordinates of all nodes
!                        1: coordinate in global x-direction
!                        2: coordinate in global y-direction
!                        3: coordinate in global z-direction
!     lakonl             element label
!     konl(1..20)        nodes belonging to the element
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     nmpc               number of MPC's
!     ikmpc(1..nmpc)     ordered global degrees of freedom of the MPC's
!                        the global degree of freedom is
!                        8*(node-1)+direction of the dependent term of
!                        the MPC (direction = 0: temperature;
!                        1-3: displacements; 4: static pressure;
!                        5-7: rotations)
!     ilmpc(1..nmpc)     ilmpc(i) is the MPC number corresponding
!                        to the reference number in ikmpc(i)   
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!     sti(i,j,k)         actual Cauchy stress component i at integration
!                        point j in element k. The components are
!                        in the order xx,yy,zz,xy,xz,yz
!     xstateini(i,j,k)   value of the state variable i at integration
!                        point j in element k at the beginning of the
!                        present increment
!     xstateini(i,j,k)   value of the state variable i at integration
!                        point j in element k at the end of the
!                        present increment
!     nstate_            number of state variables
!     dtime              time length of the increment
!
!
!     OUTPUT:
!
!     flux(1)            magnitude of the flux
!     flux(2)            not used; please do NOT assign any value
!     iscale             determines whether the flux has to be
!                        scaled for increments smaller than the 
!                        step time in static calculations
!                        0: no scaling
!                        1: scaling (default)
!           
      implicit none
!
      character*8 lakonl
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,konl(20),ipompc(*),nstate_,i,
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),node,idof,id,iscale,mi(*)
!
      real*8 flux(2),time(2),coords(3),sol,temp,press,vold(0:mi(2),*),
     &  area,co(3,*),coefmpc(*),sti(6,mi(1),*),xstate(nstate_,mi(1),*),
     &  xstateini(nstate_,mi(1),*),dtime,Start

      real*8 u(3), v(3), prod(3), norm(3), normdef(3), len, freq, so
!
      intent(in) sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi,sti,
     &     xstateini,xstate,nstate_,dtime
!
      intent(out) flux,iscale
!
!     Start of your own code.
!     NOTE from Eric - this file goes in the CalculiX 'src' directory,
!     and thus any changes here would lead to re-compilation
!     of the CCX source code.
!
!     Please note some of the useful variables (like time, etc) at the top of
!     this source 
!
!     To invoke this user defined flux in the input deck,
!     you'd need to set:
!
!     *Dflux
!      Surface_Name, S2NUx
!
!      Compute normals by (v2-v1)x(v3-v1),
!      where v2-v1 = <x2-x1,y2-y1,z2-z1>
       u(1) = co(1,konl(3)) - co(1,konl(2))
       u(2) = co(2,konl(3)) - co(2,konl(2))
       u(3) = co(3,konl(3)) - co(3,konl(2))
       
       v(1) = co(1,konl(1)) - co(1,konl(2))
       v(2) = co(2,konl(1)) - co(2,konl(2))
       v(3) = co(3,konl(1)) - co(3,konl(2))
       
       prod(1) = u(2)*v(3) - u(3)*v(2)
       prod(2) = u(3)*v(1) - u(1)*v(3)
       prod(3) = u(1)*v(2) - u(2)*v(1)

       len = sqrt(prod(1)*prod(1) + prod(2)*prod(2) + prod(3)*prod(3))

       norm(1) = prod(1)/len
       norm(2) = prod(2)/len
       norm(3) = prod(3)/len
   
!      Print to screen (to test norms are correct-ish)
!       write(*,*) norm(1), norm(2), norm(3)

!      Get deflected norms (calculating element edge vectors)
       u(1) = (co(1,konl(3)) + vold(1,konl(3))) - (co(1,konl(2)) + 
     &       vold(1,konl(2)))
       u(2) = (co(2,konl(3)) + vold(2,konl(3))) - (co(2,konl(2)) + 
     &       vold(2,konl(2)))
       u(3) = (co(3,konl(3)) + vold(3,konl(3))) - (co(3,konl(2)) + 
     &       vold(3,konl(2)))

       v(1) = (co(1,konl(1)) + vold(1,konl(1))) - (co(1,konl(2)) + 
     &       vold(1,konl(2)))
       v(2) = (co(2,konl(1)) + vold(2,konl(1))) - (co(2,konl(2)) + 
     &       vold(2,konl(2)))
       v(3) = (co(3,konl(1)) + vold(3,konl(1))) - (co(3,konl(2)) + 
     &       vold(3,konl(2)))

!      Cross product
       prod(1) = u(2)*v(3) - u(3)*v(2)
       prod(2) = u(3)*v(1) - u(1)*v(3)
       prod(3) = u(1)*v(2) - u(2)*v(1)

!      Calculate the length of the current normal vector.
       len = sqrt(prod(1)*prod(1) + prod(2)*prod(2) + prod(3)*prod(3))

!      Normalize the deflected normal direction vector components.
       normdef(1) = prod(1)/len
       normdef(2) = prod(2)/len
       normdef(3) = prod(3)/len

!      Uncomment the following line to see the deflected normals of
!      each element every iteration. 1, 2, and 3 represent x,y,z dir.
!       write(*,*) normdef(1), normdef(2), normdef(3)

!      Flux applied to the surface is modified by the dot product of the
!      reference normal and the deflected normal in the "y" direction.
!      Flux amplitude is set at 1300 W/m^2 and is centered about zero.
!      The max load +1300 W/m^2 and the min load -1300 W/m^2.
       so = 1300
!      Freqency of the loading, represents the rotation of the beam in space.

!       write(*,*) "noel=", noel
!      if (kinc .EQ. 1) then
!	 open(UNIT=100, FILE="freq")
!       	 read(100,*) Freq
!         write(*,*) Freq
!         rewind 100
!         close(100)
!        endif


!      The following line should be enabled to apply cyclic loading with a Freq input in Hz,
!      which simulates a rotating beam in a static source of irradiation.
!       Freq = 1.820460 ! Takes a frequency value in Hz to be converted into rads/s.
!      Base flux according to rigid body radiation as a function of incident
!      angle in space.
!       flux(1) = so * sin(Freq * 2 * 3.14159265 * time(2)) ! *Needs a phase change for cos(x)



!      The following lines should be enabled to apply step loading at a Start input time,
!      which simulates a sudden exposure of a beam to a static source of irradiation.
       Start = 0.5 ! Dictates the total time value where the step load is applied. 
          IF (time(2) > Start) THEN
             flux(1) = so
!      If the start time has not yet been reached, no therma flux is applied.
          ELSE
             flux(1) = 0
          ENDIF



!      Accounting for two sides, 180 degrees out of phase.
!      If the element outer surface normal vector has the same sign as the
!      oncoming flux, that corresponds to that side being irradiated.
          IF (normdef(2) .GT. 0 .AND. flux(1) .GT. 0) THEN
             flux(1) = flux(1)
          ELSEIF (normdef(2) .LT. 0 .AND. flux(1) .LT. 0) THEN
             flux(1) = ABS(flux(1))
!      If the normal vector and flux don't have the same sign, then the
!      element is shadowed and recieves no flux.
          ELSE 
             flux(1) = 0
          ENDIF

!      Print the flux based only on the rotation (A*cos(wt)).
!       write(*,*) "Base flux: ", flux(1)
!      Modify the base flux depending on axial deflection of beam.
       flux(1) = flux(1) * ABS(normdef(2))
!      Modify the flux depending on radial orientation of element.
       flux(1) = flux(1) * SQRT(1 - (normdef(3)*normdef(3)))
!      Print the adjusted flux.
!       write(*,*) "Flux adjusted for deflection: ", flux(1)

       return
       end

