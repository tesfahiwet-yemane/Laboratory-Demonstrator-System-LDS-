
!*************************************************************************************
!    Laboratory Demonstrator System (LDS)                                            *
!    Copyright (c) 2021 Geophysics Research Laboratory (GRL)                         *
!                                                                                    *
!    LDS is free software: you can redistribute it and/or modify                     *
!    it under the terms of the GNU General Public License as published by            *
!    the Free Software Foundation, either version 3 of the License, or               *
!    (at your option) any later version.                                             *
!                                                                                    *
!    LDS is distributed in the hope that it will be useful,                          *
!    but WITHOUT ANY WARRANTY; without even the implied warranty of                  *
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                   *
!    GNU General Public License for more details.                                    *
!                                                                                    *
!    You should have received a copy of the GNU General Public License               *
!    along with LDS.  If not, see <https://www.gnu.org/licenses/>.                   *
!*************************************************************************************

!*************************************************************************************
!        FORTRAN 95 program containing subroutines and other sub-programs            *
!                   for the Laboratory Demonstrator System (LDS)                     *
!  The program computes current distributions and other items for a given electrode  *
!             spacing on the basis of a programming model of the LDS                 *
!                                                                                    *
!             Geophysics Research Laboratory (GRL)                                   *
!                 Department of Earth Sciences                                       *
!             Eritrea Institute of Technology (EIT)                                  *
!                         Eritrea                                                    *
!*************************************************************************************

    module color
       implicit none
!******************************************************************************
!These are default parameters which control line width, colors, and terminal
!******************************************************************************
    character(len=3),parameter     :: default_linewidth='1'
    character(len=100),parameter   :: default_color1='blue'
    character(len=100),parameter   :: default_color2='dark-green'
    character(len=100),parameter   :: default_color3='orange-red'
    character(len=100),parameter   :: default_color4='dark-salmon'
    character(len=100),parameter   :: default_terminal='wxt'
    character(len=100),parameter   :: default_palette='CMY'
!******************************************************************************
    interface contour
        module procedure inputP
        module procedure Kfactor
        module procedure Matrix
        module procedure gaussj
        module procedure currentV
        module procedure Intensity
        module procedure car
        module procedure twoD
    end interface
!******************************************************************************
    contains
!******************************************************************************
    function my_date_and_time() result(f_result)
!******************************************************************************
!This function creates a string with current date and time
!It is a default method to name output files, if needed
!******************************************************************************
    implicit none
    character(len=8)   :: date
    character(len=10)  :: time
    character(len=33)  :: f_result
!******************************************************************************
    call date_and_time(date,time)
    f_result='date_'//date(7:8)//'-'//date(5:6)//'-'//date(1:4)//'_time_'//time(1:2)//':'//time(3:4)//':'//time(5:10)
!******************************************************************************
    end function my_date_and_time
!******************************************************************************
    function output_terminal(terminal) result(f_result)
!******************************************************************************
    implicit none
    character(len=*),intent(in) :: terminal
    integer, parameter          :: Nc=35
    character(len=Nc)           :: f_result
!******************************************************************************
    select case(terminal)
        case('ps')
            f_result='postscript landscape color'
            case default
                f_result=terminal
    end select
    end function output_terminal
!******************************************************************************

!******************************************************************************
!Subroutines for the main program LDS
!******************************************************************************

    subroutine InputP (IL, IR, VL, VR, NNH, NNV, d, V1, CI, disp, DEV)
    implicit none
!******************************************************************************
!This subroutine inputs the parameters that are needed by the program. There are two files that are found in the same directory
!as the program that are read. These are the files 'config' which contains the parameters and the 'profile' file that contains
!the weight of the profiles to be considered.
!******************************************************************************
!Declaration of Variables
!******************************************************************************
    integer  :: IR, IL, VR, VL, NNH, NNV, DEV, i, CI(30), disp, w
    real     :: d, j
    real(kind = 8) V1
!******************************************************************************
!Reading the input parameters from a file 'config'
!******************************************************************************
    open(unit = 11, file = 'config', status = 'old')
    read(11, *) IL, IR, VL, VR, NNH, NNV, d, V1
!13 format (6(2x,I5), 2(2x,F4.2))
    write(*,*)'******************************************************************************'
    write(*,*)'          The Laboratory Demonstrator System (LDS) '
    write(*,*)'******************************************************************************'
    write(*,*)'IL  :',IL
    write(*,*)'IR  :',IR
    write(*,*)'VL  :',VL
    write(*,*)'VR  :',VR
    write(*,*)'NNH :',NNH
    write(*,*)'NNV :',NNV
    write(*,*)
    close(11)
!******************************************************************************
    write(*,*)'******************************************************************************'
9   write(*,*) 'If Wenner config is desired, input 1'
    write(*,*)'******************************************************************************'
    read (*,*) i
    write(*,*)'******************************************************************************'
    if (i.eq.1) goto 5
!******************************************************************************
!Make runtime modifications to the initial input for Schulmberger config
!******************************************************************************
    write(*,*)'******************************************************************************'
    write(*,*)'Input deviation from the initial setting'
    write(*,*)'******************************************************************************'
    read(*,*) DEV
    IL = IL - DEV
    IR = IR + DEV
    write(*,*) '******************************************************************************'
    write(*,*) 'input Horizontal displacement'
    read(*,*) disp
    IL = IL + disp
    IR = IR + disp
    VL = VL + disp
    VR = VR + disp
    write(*,*)'******************************************************************************'
    write(*,*)'IL  :',IL
    write(*,*)'IR  :',IR
    write(*,*)'VL  :',VL
    write(*,*)'VR  :',VR
    write(*,*)'******************************************************************************'
    write(*,*)
    write(*,*) 'If correction is desired enter 1 and input correct values, or another number to continue'
    read(*,*) i
    if(i.eq.1) goto 9
    goto 6
!******************************************************************************
!Determine spacing of the 4 electrodes for the Wenner method
    write(*,*)'******************************************************************************'
5 write(*,*) 'Input deviation from the initial setting'
    write(*,*)'******************************************************************************'
    read(*,*) DEV
    j = DEV/2.0
    i = int(j)
    j = j - i
    if(j.gt.0.or.DEV.gt.16) then
        write(*,*) 'Entered value must be a positive even number less than 17; please try again'
        goto 5
    end if
    write(*,*) '******************************************************************************'
    write(*,*) 'input Horizontal displacement'
    read(*,*) disp
    VL = VL - DEV/2
    VR = VR + DEV/2
    IL = 2*VL - VR
    IR = 2*VR -VL
    VL = VL + disp
    VR = VR + disp
    IL = IL + disp
    IR = IR + disp
    write(*,*)'******************************************************************************'
    write(*,*)'IL  :',IL
    write(*,*)'IR  :',IR
    write(*,*)'VL  :',VL
    write(*,*)'VR  :',VR
    write(*,*)'******************************************************************************'
    write(*,*)
    write(*,*) 'If correction is desired enter 1 and input correct values, or another number to continue'
    read(*,*) i
    if(i.eq.1) goto 5
    write(*,*)'Press return to continue'
    read(*,*)
    goto 6
    write(*,*)'******************************************************************************'
    write(*,*)'IL  :',IL
    write(*,*)'IR  :',IR
    write(*,*)'VL  :',VL
    write(*,*)'VR  :',VR
    write(*,*)'******************************************************************************'
    write(*,*)
    write(*,*) 'If correction is desired enter 1 and input correct values, or another number to continue'
    read(*,*) i
    if(i.eq.1) goto 9
    write(*,*)'Press return to continue'
    read(*,*)
!******************************************************************************
!Input the profile weights
6   open(unit=15, file='profile',status = 'old')
    read(15,*) CI
    close(15)
!******************************************************************************
    return
!******************************************************************************
end subroutine InputP
!******************************************************************************

subroutine Kfactor(IL, IR, VL, VR, Kfact, d)
    implicit none
!******************************************************************************
! This subroutine computes the geometric factor, K for a given electrode arrangement
! This is done only for comparison and general teaching purposes and does not influence
! the computational aspect of the system
!******************************************************************************
!IL >>> Position of the left current electrode
!IR >>> Position of the right current electrode
!VL >>> Position of the left voltage electrode
!VR >>> Position of the right voltage electrode
!PA = rA, PB = rB
!******************************************************************************
!Declaration of Variables
!******************************************************************************
    integer                 :: IL, IR, VL, VR
    real(kind=8)            :: Kfact, Temp
    real                    :: d, PA, PB, RA, RB
    real(kind=8), parameter :: pi = 4*atan(1.0)
!******************************************************************************
    PA = d*(VL - IL)
    PB = d*(IR - VL)
    RA = d*(VR - IL)
    RB = d*(IR - VR)
!Calculating the denominator of the Kfactor
    Temp = ((1.0/PA - 1.0/PB)-(1.0/RA - 1.0/RB))
!The final Kfactor
    Kfact = (2.0*pi)/Temp
    write(*,*)'******************************************************************************'
    write(*,*)'Kfactor   :',Kfact
    write(*,*)'******************************************************************************'
    return
end subroutine Kfactor
!******************************************************************************
    subroutine Matrix(IL, IR, NNH, NNV, N, M, CI, cmp)
    implicit none
!******************************************************************************
!This subroutine generates the matrix that is formed by writing and re-arranging the loop equations of the demonstrator model
!circuit. The resulting equation is of the form MI=V, where M is the matrix while I and V are the transposes of the loop
!currents and source vectors respectively.
!Modification to introduce two high resistance spots. The values are 1000 times greater than the other resistors. Done 20/05/2018
!Modification to accommodate truncation effect introduced on 01/02/2018
!Only the subroutine Matrix has been affected
!******************************************************************************
!Declaration of variables
!******************************************************************************
    integer       :: IL, IR, NNH, NNV, TP, R, i, j, CI(30)
    integer       :: CU, RL, CL, RU, N
    real(kind=8)  :: M(1422,1422)
    real          :: AV, CMP !CMP is the compensation factor for correcting the truncation effect
!******************************************************************************
    write(*,*)'******************************************************************************'
15  write(*,*)'Input the compensation factor: It should be 0.75 < CMP < 1  '
    write(*,*)'******************************************************************************'
    read(*,*) CMP
    write(*,*)'******************************************************************************'
    if(CMP.LT.0.75) goto 15
    if(CMP.GT.1) goto 15
!******************************************************************************
!Set the value of M to 0 (clear M). The matrix has a dimension of N x N, where N = (NNH-1)*(NNV-1) + 1
    N = (NNH-1)*(NNV-1)+ 1
    do i = 1, N
        do j = 1, N
            M(i,j)= 0
        end do
    end do
!******************************************************************************
!Determine the first row of the matrix
    TP = IR - IL + 1
!This is essentially the number of segments traversed + 1 (1 is added to cover the connecting segments
    M(1,1) = TP*CI(1)
    do i = IL, IR - 1
        M(1,i+1) = -1*CI(1)
    end do
!******************************************************************************
!Determine the rows of the matrix produced by the loops in the first row of the model
!******************************************************************************
    AV = (CI(1) + CI(2))/2
    if (IL.gt.1)then
        do i = 1, IL - 1
            R  = i + 1
            CU = R + 1
            RU = NNH + i
            M(R,R) = 4*AV
            M(R,i) = -1*CI(1)
            if (i.eq.1) M(R,i) = 0
            M(R,RU) = -1*CI(2)
            M(R,CU) = -1*AV
        end do
    end if
    do i = IL, IR-1
        R  = i + 1
        CU = R + 1
        RU = NNH + i
        M(R,R) = 4*AV
        M(R,i) = -1*CI(1)
        M(R,CU)= -1*AV
        M(R,RU)= -1*CI(2)
        M(R,1) = -1*CI(1)
    end do
    do i = IR, NNH-1
        R  = i + 1
        CU = R + 1
        RU = NNH + i
        M(R,R) = 4*AV
        M(R,i) = -1*CI(1)
        M(R,CU)= -1*AV
        M(R,RU)= -1*CI(2)
        if(i.eq.NNH-1) M(R,CU)=0
    end do
!******************************************************************************
!Determine the elements of the rows of the matrix corresponding to the remaining rows of the model
!******************************************************************************
    do j = 2, NNV - 1
        AV = (CI(j)+ CI(j+1))/2
        do i = 1, NNH - 1
            R  = NNH + (j - 2)*(NNH - 1) + i
            CU = R + 1
            CL = R - 1
            RU = R + NNH - 1
            RL = R - NNH + 1
            M(R,R) = 4*AV
            M(R,CL)= -1*AV
            !write (*,*) CU, CL, RU, RL, i, j
            if(i.eq.1) M(R,CL) = 0
            if(RU.gt.N) goto 20
            M(R,RU) = -1*CI(j+1)
            20 M(R,RL) = -1*AV
            if(CU.gt.N) goto 25
            M(R,CU) = -1*AV
            !write(*,*)'RU = ', RU, '    CU = ', CU
            if(i.eq.NNH-1) M(R,CU) = 0
      25  end do
    end do
!35  write(*,*) 'input the row and columns separated by a comma'
!    read(*,*) i, j
!    write(*,*) M(i,j)
!    write(*,*)
!    goto 35
!******************************************************************************
!!Input high resistance locations
!    do j = 10, 25
!        AV = 1000*(CI(j) + CI(j+1))/2
!        do i = 10, 25
!            R = NNH + (j-2)*(NNH-1) + i
!            CU = R + 1
!            CL = R - 1
!            RU = R + NNH - 1
!            RL = R - NNH + 1
!            M(R,R) = 4*AV
!            M(R,CL) = -1*AV
!!write (*,*) cu, cl, ru, rl, i, j
!            if(i.eq.1) M(R,CL) = 0
!            if(RU.gt.N) goto 18
!            M(R,RU) = -1000*CI(j+1)
!            18 M(R,RL) = -1000*AV
!            if(CU.gt.N) goto 23
!            M(R,CU) = -1*AV
!        23 end do
!    end do
!!******************************************************************************
!!put a second high resistance spot
!    do j = 10, 25
!        AV = 1000*(CI(j) + CI(j+1))/2
!        do i = 30, 45
!            R = NNH + (j-2)*(NNH-1) + i
!            CU = R + 1
!            CL = R - 1
!            RU = R + NNH - 1
!            RL = R - NNH + 1
!            M(R,R) = 4*AV
!            M(R,CL) = -1*AV
!!write (*,*) cu, cl, ru, rl, i, j
!            if(i.eq.1) M(R,CL) = 0
!            if(RU.gt.N) goto 19
!            M(R,RU) = -1000*CI(j+1)
!            19 M(R,RL) = -1000*AV
!            if(CU.gt.N) goto 24
!            M(R,CU) = -1*AV
!        24 end do
!    end do
!******************************************************************************
! Make corrections for the truncation effect
    M(2,2) = M(2,2)*CMP
    do i = 1, NNV-2
        j = i*(NNH-1) + 1
        M(j,j) = M(j,j)*CMP
        M(j+1,j+1) = M(j+1,j+1)*CMP
!   write(*,*) ' i = ' ,i, ' j = ' ,j, ' M(j,j) = ' , M(j,j), ' M(j+2, j+2) ' M(j+2,j+2) = ', ' N = ',N
    end do
    j = j + 2 !This part takes care of the last row of the programming model
    do i = j, N
        M(i,i) = M(i,i)*CMP
    end do
!   write(*,*) M(1375,1375)
!   10 write(*,*) 'Input row number'
!   read(*,*) i
!   do j = 1, N
!   write(*,*) M(i,j)
!   end do
!   goto 10
!   write(*,*) 'matrix complete'
!******************************************************************************
end subroutine Matrix
!******************************************************************************

subroutine gaussj(a, b, n)
    implicit none
!******************************************************************************
!This subroutine inverts the matrix. This subroutine is adapted from the book
!'Numerical recipes: The art of scientific computing, third ed.' Cambridge
!University Press; by Press et al.
!******************************************************************************
!Declaration of variables
!******************************************************************************
    integer              :: q, n, NMAX
    real(kind=8)         :: b(1422,1422), a(1422,1422)
    parameter (NMAX = 1422)
    integer              :: i, icol, irow, j, k, l, ll, indxc(NMAX), indxr(NMAX)
    integer              :: ipiv(NMAX)
    real(kind=8)         :: dum, pivinv
    real(kind=8)         :: big
!******************************************************************************
    q = n
    !a = M
    !b = MI
    !write(*,*) n,m
    !write(*,*)
!******************************************************************************
    do j = 1, n
        ipiv(j) = 0
    end do
    do i = 1,n
        big = 0
        do j = 1,n
            if(ipiv(j).ne.1)then
                do k = 1,n
                if (ipiv(k).eq.0)then
                if(abs(a(j,k)).ge.big)then
                    big = abs(a(j,k))
                    irow = j
                    icol = k
                end if
                elseif(ipiv(k).gt.1)then
                    continue !'singular matrix 1'
            end if
        end do
    end if
end do

!******************************************************************************
    ipiv(icol) = ipiv(icol)+ 1
    if(irow.ne.icol)then
        do l = 1, n
            dum = a(irow,l)
            a(irow,l) = a(icol,l)
            a(icol,l) = dum
        end do
        !write(*,*)'1'
        do l = 1, q
            dum = b(irow,l)
            b(irow,l) = b(icol,l)
            b(icol,l) = dum
        end do
    end if
!******************************************************************************
    !write(*,*) 'tt'
    indxr(i) = irow
    indxc(i) = icol
    !write(*,*) a(icol,icol), icol, irow
    if(a(icol,icol).eq.0) continue !pause 'Singular matrix 2'
    pivinv = 1./a(icol,icol)
    a(icol,icol) = 1
    do l = 1, n
        a(icol,l) = a(icol,l)*pivinv
    end do
!******************************************************************************
    do l = 1,q
        b(icol,l) = b(icol,l)*pivinv
    end do
    do ll = 1,n
        if(ll.ne.icol)then
            dum = a(ll,icol)
            a(ll,icol) = 0
            do l = 1, n
                a(ll,l) = a(ll,l) - a(icol,l)*dum
            end do
            do l = 1,q
                b(ll,l) = b(ll,l) - b(icol,l)*dum
            end do
        end if
        end do
    end do
!*************************
!*****************************************************
    do l =  n, 1, -1
        if(indxr(l).ne.indxc(l))then
            do k = 1,n
                dum = a(k,indxr(l))
                a(k,indxr(l))= a(k,indxc(l))
                a(k,indxc(l))= dum
            end do
        end if
    end do
!******************************************************************************
    !10 write(*,*) 'Input column number'
    !read(*,*) i
    !do j = 1, N
    ! write(*,*) a(j,i)
    !end do
    !go to 10
end subroutine gaussj
!******************************************************************************
subroutine currentV(NNH, NNV, M, V1, IC, I0)
    implicit none
!******************************************************************************
!This subroutine computes the values of all currents flowing in all the branches of the LDS programming model. It gets the
!following inputs: The number of nodes in the horizontal direction of the model NNH, the number of nodes in the vertical direction
!of the model NNV, the inverted matrix M, the value of the battery voltage in the model V1, It gets as an output the current
!vector, which has naming in the form of matrix I(i,j)
!******************************************************************************
!Declaration of variables
!******************************************************************************
    integer                     :: NNH, NNV, K, N
    integer                     :: i, j
    real(kind=8)                :: M(1422,1422), IC(1422,1422), I0, V1
!******************************************************************************
!Get the value of the currents
!******************************************************************************
    !open(unit=13, file='current.txt', status='replace')
!******************************************************************************
!Clear the matrix IC
    N = (NNV - 1)*(NNH - 1) + 1
    do j = 1, N
        do i = 1, N
            IC(j,i) = 0
        end do
    end do
!******************************************************************************
    open(unit= 28, file='currents', status='replace')
    I0 = M(1,1)*V1
    K  = 1
    do j = 1, NNV - 1
        do i = 1, NNH - 1
            K = K + 1
            IC(j,i) = M(1,K)*V1
            write(28,*) i, j, IC(j,i)
        end do
    end do
    close(28)
!******************************************************************************
    write (*,*)
    write(*,*)'******************************************************************************'
    write(*,*) 'I0      :', I0
    write(*,*)'******************************************************************************'
    !do i = 1, K
    !IC(2,i) = M(1,i)*V1
    !write(*,*) IC(2,i)
    !end do
    !return
    return

end subroutine currentV
!******************************************************************************
subroutine Intensity(IC, NNH, NNV, MT, IL, IR, I0)
!First written on 31 December 2017
!This subroutine obtains the current values flowing through every branch of the model.
!The result is stored in a matrix called MT, and this will serve as an input to a graphic display subroutine
!******************************************************************************
!Declaration of variables
!******************************************************************************
    integer      :: IL, IR, NNH, NNV, N, DEV, i, j, k
    real(kind=8) :: IC(1422,1422), MT(1422,1422), I0
!******************************************************************************
!Clear the matrix MT
    N = (NNH-1)*(NNV-1) + 1
    do j = 1, N
        do i = 1, N
            MT(j,i) = 0
        end do
    end do
!Determine the current values in the first row of the model
    open(unit=13, file='Intensity.txt', status='unknown')
    do i = 1, IL-1
        MT(1,i) = 0 - IC(1,i)
        write(13,*) i, 1, MT(1,i)*1000
    end do
    do i = IL, IR - 1
        MT(1,i) = -IC(1,i) + I0
        write(13,*) i, 1, MT(1,i)*1000
    end do
    do i = IR, NNH - 1
        MT(1,i) = 0 - IC(1,i)
        write(13,*) i, 1, MT(1,i)*1000
    end do
!Determine the currents of the values of the other rows
    do j = 2, NNV
        k = j-1
        do i = 1, NNH-1
            MT(j,i)= -IC(j,i) + IC(k,i)
            write(13,*) i, j, MT(j,i)*1000
        end do
    end do
    close(13)
    DEV = NNH/2
    do i = 1, NNV
        write(*,*) i, '            ', MT(i,Dev)
    end do

!******************************************************************************
end subroutine Intensity
!******************************************************************************

subroutine car(VL, VR, IC, I0, Kfact, AR, NNH, NNV, IR, IL, CI, disp, dater, DEV, V1, MT, cmp)
    implicit none
!******************************************************************************
!This subroutine computes the first factor of the apparent resistivity, that is the term delta V(dV) over I0.
!The values received by this subroutine are: The current I0, The position of the potential electrodes VL and VR, the value
!resistive arms Kfact, when this is available. The result obtained here is intended for comparison purposes and does't
! affect system operation
!******************************************************************************
    integer      :: VL, VR, i, C, DEV, NNH, NNV, K1, IR, IL, DE, CI(30), disp, j, DEVV, ans
    real         :: cmp
    real(kind=8) :: I0, Kfact, dV, IC(1422,1422), AR, CR(30), App(30), Total, CI0(30),V1, CR0(30), MT(1422,1422)
    real(kind=8) :: TUP(30), TDW(30), DPper, DSper(30), temp, tt(30), di(30), rr, rrt, sums !temporary for CR
    character*10 :: dater*12, report*16, temporary*10
    character*1  :: yn
!******************************************************************************
!The current through the voltmeter in the programming model is assumed to be zero (0)
!compute dV
!******************************************************************************
    dV = 0
    do i = VL, VR - 1
        dV = dV + (I0 - IC(1,i))*CI(1) !Determine the reading of the voltmeter
    end do
    AR = (Kfact*dV)/I0
    write(*,*)'******************************************************************************'
! Modification of the program to write the apparent resistivity in a file on 16/12/2018
    open(unit=17, file='Apparent.txt', status='unknown')
    write(17, *) AR
    write(*,*)'******************************************************************************'
    write(*,*) 'Potential difference:', dV
    write(*,*)'******************************************************************************'
!Print the values of the currents in the middle of the model, or finish
    write(*,*)
    write(*,*)'For printing currents press 1, or from 2 - 9 to finish'
    write(*,*)'******************************************************************************'
    read (*,*) C
    if(C.ne.1.)goto 10
!Compute the values of the currents
    DEVV   = (NNH/2) + disp
    write(*,*) DEVV, disp
    CR(1) = I0 - IC(1,DEVV)
    Total = CR(1)
    do i = 1, NNV
        K1 = i + 1
        CR(i+1) = IC(i,DEVV) - IC(K1,DEVV)
        Total = Total + CR(i+1)
    end do
    write(*,*) 'Input the desired output mode: 0 for normal mode, 1 for defects mode    '
    read(*,*) ans
    if(ans.eq.1) then
        call defects(IL,IR,VL,VR,V1,dater,I0,cmp,DEV,disp,CR,CR0,NNH,NNV,CI,MT)
        goto 50
    endif
    do j = 1, NNV     ! This computes the percentage of the total current flowing
        TUP(j) = 0    ! above (TUP) and below (TDW) a certain layer
        do i = 1, j
            TUP(j) = TUP(j) + CR(i)
        end do
        TUP(j) = TUP(j)*1000 ! Values in ma
        TDW(j) = 0
        do i = j+1, NNV
            TDW(j) = TDW(j) + CR(i)
        end do
        TDW(j) = TDW(j)*1000 !Values in ma
    end do
     if(disp.eq.0) then
          open(unit = 24, file = 'DEV0', status = 'new')
          do i = 1, NNV
              write(24,*) CR(i)
          end do
          close(24)
     else
          open(unit = 25, file = 'DEV0', status = 'old')
          do i = 1, NNV
              read(25,*) CR0(i)
          end do
          close(25)
     end if
     write(*,*) 'input the name of the file for the report'
     read(*,*) report
     open(unit = 19, file = report, status = 'new')
     write(19,*) '                        ELECTRICAL RESISTIVITY LDS MODEL OUTPUTS'
     write(19,12) DEV, disp, cmp, dater
     rr = I0*1000
     write(19,11) IL, IR, VL, VR, rr, V1
     write(19,*)
     if(disp.eq.0) write(19,*) ' Layer  CI(Ohms)      CR(ma)        Tup(ma)       Tdw(ma)      up%'
     if(disp.ne.0) write(19,*) ' Layer  CI(Ohms)       CR(ma)     Tup(ma)     Tdw(ma)     up%         di%'
     !write(19,*)
     rrt = 0
     j = 0
17 sums = 0
     do i = 1, NNV
        di(i) = 1000000*(CR(i) - CR(i+1))
        if(i.eq.30) rrt = CR(i)
        if(i.eq.29.or.i.eq.30) then
            sums = sums + 0
        else
            sums = sums + di(i) - rrt
        end if
        end do
        if(j.eq.0) then
            j = 1
            goto 17
        end if
        do i = 1, NNV
            DPper = (TUP(i)/I0)/10		!values in percentage
        tt(i) = CR(i)*1000
        if(disp.ne.0) DSper(i) = 100*(CR(i) -CR0(i))/CR0(i)
        write(*,*) CR(i), CR0(i), DSper(i)
        di(i) = (di(i)/sums)*100
        !tt = CR(i)*1000
        !if(disp.ne.0) DSper = 100*(CR(i) - CR0(i))/CR0(i)
        !write(*,*) CR(i), CR0(i), DSper
        if(disp.eq.0) write(19,13) i, CI(i), tt(i), TUP(i), TDW(i), DPper
        if(disp.ne.0) write(19,14) i, CI(i), tt(i), TUP(i), TDW(i), DPper,      di(i)
     end do
     close(19)
     write(*,*) sums
 11  format(2x, 'IL = ', i2, 3x, 'IR = ', i2, 3x, 'VL = ', i2, 3x, 'VR = ', i2, 3x, 'I0 (ma) = ', F7.5, 3x, 'Bat(vlt) = ', F7.2)
 12  format(6x, 'Deviation: ', i2, 6x, 'Displacement: ', i2, 6x, 'Cmp fact: ', F5.3, 6x, 'Date: ', A12)
 13  format(2x, i2,6x, i5, 6x, F8.4, 6x, F8.4, 6x, F8.4, 6x, F7.3)
 14  format(2x, i2,6x, i5, 6x, F8.4, 4x, F8.4, 4x, F8.4, 4x, F8.4, 4x, F8.4)
50   write(*,*)
    write(*,*)'******************************************************************************'
    write(*,*)'******************************************************************************'
    write(*,*)'I0      : ',I0
    write(*,*)'Total   : ',Total
    write(*,*)'IL      : ',IL
    write(*,*)'IR      : ',IR
    write(*,*)'VL      : ',VL
    write(*,*)'VR      : ',VR
    write(*,*)'Kfactor : ',Kfact
    write(*,*)'______________________________________________________________________________'
!******************************************************************************
!Print the values of the currents
!******************************************************************************
    write(*,*) '       LAYER NUMBER                              VALUE OF CURRENT'
    write(*,*)'______________________________________________________________________________'
    do i = 1, NNV
        write(*,*)'  ',i, '                            ',CR(i)
    end do
    write(*,*)'______________________________________________________________________________'
!******************************************************************************
!    open(unit=10, file='test.txt', status='replace')
!     do j = 1, NNV-1
!        do i =1, NNH - 1
!    write(10, *) i,j,CR(i)
!    end do
!end do
!******************************************************************************
!close(unit=13)
 75   write(*,*) 'Do you want to write current values? (y/n)'
      read(*,*) yn
      if(yn.eq.'n') go to 10
      write(*,*)'Input column number'
      read(*,*) j
      write(*,*) 'Column number  =  ', j
      write(*,*) '  Layer         current value'
      sums = 0
      do i = 1, NNV
         write(*,77) i, MT(i,j)
         sums = sums + MT(i,j)
      end do
      write(*,*)
      write(*,*) ' Sum = ', sums
      goto 75
 77   format(2x, i2, 3x, F16.14)
10    return
end subroutine car
!*******************************************************************************************************
     subroutine defects(IL,IR,VL,VR,V1,dater,I0,cmp,DEV,disp,CR,CR0,NNH,NNV,CI,MT)
     implicit none
!*******************************************************************************************************
!This subroutine produces the output option 1. This option is intended to deal with the effects of the truncated
!nature of the model. It produces results that are needed to quantify the effects.
!*******************************************************************************************************
    integer      :: VL, VR, DEV, NNH, NNV, IR, IL, CI(30), disp, i, j
    real         :: cmp
    real(kind=8) :: I0, CR(30), App(30), CI0(30),V1, CR0(30), MT(1422,1422), rr   !rr is I0*1000 so that it is in ma
    real(kind=8) :: DSper(30), temp,tt(30), di(30), L(30), R(30)
    character    :: dater*12, report*30
!*******************************************************************************************************
!Save into file the current profile if DEV = 0 and ds = 0
!*******************************************************************************************************
     if(disp.eq.0.and.DEV.eq.0) then
          open(unit = 24, file = 'DEV0', status = 'replace')  !Stored in file DEV0 for reference purposes
          do i = 1, NNV
              write(24,*) CR(i)
          end do
          close(24)
     else
          open(unit = 25, file = 'DEV0', status = 'old')  ! Read to be used for comparison purposes
          do i = 1, NNV
              read(25,*) CR0(i)
          end do
          close(25)
     end if
     write(*,*) 'input the name of the file for the report'
     read(*,*) report
     open(unit = 19, file = report, status = 'new')
     write(19,*) '                    ELECTRICAL RESISTIVITY LDS MODEL OUTPUTS'
     write(19,12) DEV, disp, cmp, dater
     rr = I0*1000
     write(19,11) IL, IR, VL, VR, rr, V1
     write(19,*)
     if(disp.eq.0) write(19,*) ' Layer  CI(Ohms)      L(ma)       R(ma)        CR(ma) ' ! write title for the variables
     if(disp.ne.0) write(19,*) ' Layer  CI(Ohms)      L(ma)       R(ma)       CR0(ma)     CR(ma)      Ds%'

      do i = 1, NNV
        tt(i) = CR(i)*1000
        if(disp.ne.0) DSper(i) = 100*(CR(i) - CR0(i))/CR0(i)
        write(*,*) CR(i), CR0(i), DSper(i)
        L(i) = 0
        do j = 1, VL - 1
           L(i) = L(i) + MT(i,j)      !Sum of the left side segment currents for row i, excluding the middle column
        end do
        L(i) = L(i)*1000   ! L is given in ma
        R(i) = 0
        do j = VR, NNH - 1
           R(i) = R(i) + MT(i,j)      !Sum of the right side segment currents for row i, excluding the middle column
        end do
        R(i) = R(i)*1000    ! R is given in ma
        if(disp.eq.0) write(19,13) i, CI(i), L(i), R(i), tt(i)
        if(disp.ne.0) write(19,14) i, CI(i), L(i), R(i), CR0(i)*1000, tt(i), DSper(i)
     end do
     close(19)

 11  format(2x, 'IL = ', i2, 3x, 'IR = ', i2, 3x, 'VL = ', i2, 3x, 'VR = ', i2, 3x, 'I0 (ma) = ', F7.4, 3x, 'Bat(vlt) = ', F7.2)
 12  format(2x, 'Deviation: ', i2, 6x, 'Displacement: ', i2, 6x, 'Cmp fact: ', F5.3, 6x, 'Date: ', A12)
 13  format(2x, i2,6x, i5, 6x, F8.4, 4x, F8.4, 5x, F8.4, 6x, F7.3)
 14  format(2x, i2,6x, i5, 6x, F8.4, 4x, F8.4, 4x, F8.4, 4x, F8.4, 4x, F8.4, 3x, F8.3, 4x, F8.4)
!
!           stop and return
!
     return
     end subroutine defects

!******************************************************************************
!******************************************************************************
subroutine twoD(NNV,NNH,CR, pause,palette,terminal,filename,pm3d,contour,persist,input)
    implicit none
!******************************************************************************
!This subroutine computes the color display of the current distribution
!******************************************************************************
    real, optional              :: pause
    character(len=*),optional   :: palette,terminal,filename,pm3d,contour,persist,input
    integer                     :: nx,ny,nz, NNH, NNV
    integer                     :: r,s
    integer                     :: ierror, ios, file_unit
    character(len=100)          :: data_file_name, command_file_name, my_pause, my_persist
    real(kind=8), intent(in)    :: CR(30)!x,y,z
!******************************************************************************
    nx=NNH-1
    ny=NNV-1
    nz=nx*ny
!******************************************************************************
!******************************************************************************
   if(present(input))then
        data_file_name='data_file_'//input//'.txt'
        command_file_name='command_file_'//input//'.txt'
        else
            data_file_name='data_file.txt'
            command_file_name='command_file.txt'
    end if
!******************************************************************************
    ierror=0
    call get_unit(file_unit)
    if(file_unit==0)then
        ierror=1
        print *,'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    end if
    open(unit=file_unit, file=data_file_name, status='replace',iostat=ios)
    if(ios/=0)then
        ierror=2
        print *, 'write_vector_data - fatal error! Could not open the terminal data file.'
        stop
    end if
!******************************************************************************
!Here we write the date to the data_file - the gnuplot will read this data later
!******************************************************************************
    do s = 1,nx
        do r = 1,ny
            write(file_unit,'(I5,I5,I5)') s,r,r*s
        end do
        write(file_unit,'(a)')
    end do
!******************************************************************************
    close(unit=file_unit)
!******************************************************************************
    ierror=0
    call get_unit(file_unit)
    if(file_unit==0)then
        ierror=1
        print *, 'write_vector_date - fatal error! Could not get a free FORTRAN unit.'
        stop
    end if
    open(unit=file_unit, file=command_file_name, status='replace',iostat=ios)
    if(ios/=0)then
        ierror=2
        print *, 'write_vector_data - fatal error! Could not open the terminal command file.'
        stop
    end if
!******************************************************************************
!Here we write the commands to the commands file which gnuplot will execute
!******************************************************************************
    my_persist='persist '
if (present(persist).and.(persist=='no')) my_persist=' '
if (present(terminal)) then
!    write ( file_unit, '(a)' ) 'set terminal '// trim(output_terminal(terminal))
!	if (present(filename)) then
!		write ( file_unit, '(a)' ) 'set output "'// trim(filename) //'"'
!	else
!		write ( file_unit, '(a)' ) 'set output "'//my_date_and_time()//'"'
!	end if
else
    write ( file_unit, '(a)' ) 'set terminal ' // trim(default_terminal) // ' ' &
			& //trim(my_persist) // ' title  "Gnuplot"'
end if
!******************************************************************************
    write ( file_unit, '(a)' ) 'set nokey'
if (present(palette)) then
if ((trim(palette).ne.'RGB').and.(trim(palette).ne.'HSV').and.(trim(palette).ne.'CMY').and.&
		& (trim(palette).ne.'YIQ').and.(trim(palette).ne.'XYZ')) then
    write ( file_unit, '(a)' ) 'set palette '// trim(palette)
else
    write ( file_unit, '(a)' ) 'set palette model '// trim(palette)
end if
else
    write ( file_unit, '(a)' ) 'set palette model '// trim(default_palette)
end if
!******************************************************************************
if (present(pm3d)) then
    write ( file_unit, '(a)' ) 'set '// pm3d
else
    write ( file_unit, '(a)' ) 'set pm3d map impl'
    if (present(contour)) then
        if (contour=='surface') then
            write ( file_unit, '(a)' ) 'set contour surface'
        elseif (contour=='both') then
            write ( file_unit, '(a)' ) 'set contour both'
        else
            write ( file_unit, '(a)' ) 'set contour'
        end if
    end if
end if
write ( file_unit, '(a)' ) 'set contour'
write ( file_unit, '(a)' ) 'set style increment user'
!write ( file_unit, '(a)' ) 'set style line lc rgb "black"'
!write ( file_unit, '(a)' ) 'set palette defined'
write ( file_unit, '(a)' ) 'set palette defined (0 "blue",1 "green",3 "yellow",4 "red")'
!write ( file_unit, '(a)' ) 'set palette defined (0 "blue",0.02 "green",0.1 "yellow",0.14 "red")'
!write ( file_unit, '(a)' ) 'set palette model RGB'
!write ( file_unit, '(a)' ) 'set palette defined(0 "blue",1 "dark-blue",1 "green",2 &
!&"dark-green",2 "yellow",3 "dark-yellow",4 "red",5 "dark-red")'
write ( file_unit, '(a)' ) 'set autoscale fix'
write ( file_unit, '(a)' ) 'set dgrid3d 30,50'
write ( file_unit, '(a)' ) 'set terminal wxt size 800, 500'
write ( file_unit, '(a)' ) 'set origin -0.1,-0.1'
write ( file_unit, '(a)' ) 'set size 1.2,1.2'
write ( file_unit, '(a)' ) 'set cntrparam cubicspline'
write ( file_unit, '(a)' ) 'set cntrparam points 7'
write ( file_unit, '(a)' ) 'set cntrparam order 10'
write ( file_unit, '(a)' ) 'set cntrparam levels 10'
write ( file_unit, '(a)' ) 'set cntrparam levels incr 0.002,1.2,1.5' !0.01 to 1
write ( file_unit, '(a)' ) 'set style increment userstyle'
write ( file_unit, '(a)' ) 'set style line 2 lc rgb "black"'
write ( file_unit, '(a)' ) 'set grid lw 5 ls 10'
!write ( file_unit, '(a)' ) 'set title "Laboratory Demonstrator System(LDS)" offset 0,-2.7,0 font "Arial:Bold,16"'
!write ( file_unit, '(a)' ) 'x=4'
!write ( file_unit, '(a)' ) 'y=2'
!write ( file_unit, '(a)' ) 'z=8'
!write ( file_unit, '(a)' ) 'set label "Deviation from input:" at screen 0.2,0.9 font "Arial:Bold,8"'
!write ( file_unit, '(a)' ) 'if(x<3) set label "2" at screen 0.35,0.9 font "Arial:Bold,8"; else set label "4" &
!& at screen 0.35,0.9 font "Arial:Bold,8"'
!write ( file_unit, '(a)' ) 'set label "2" at screen 0.35,0.9 font "Arial:Bold,8"'
write ( file_unit, '(a)' ) 'set xlabel "Horizontal Nodes" offset 0,-1,0 font "Arial:Bold,10"'
!write ( file_unit, '(a)' ) 'set ylabel "Number of\n Vertical Nodes" offset 0,1 rotate by 0 font "Arial:Bold,10"'
write ( file_unit, '(a)' ) 'set ylabel "Vertical Nodes" offset -6, 0 rotate parallel font "Arial:Bold,10"'
write ( file_unit, '(a)' ) 'set colorbox horiz user origin .17,.12 size .66,.07'
write ( file_unit, '(a)' ) 'set cbrange [0:0.14]'
!write ( file_unit, '(a)' ) 'set cbrange [0:0.5]'
write ( file_unit, '(a)' ) 'set cblabel "Current Intensity(mA)" offset 0,0.9,0 font "Arial:Bold,10"'
write ( file_unit, '(a)' ) 'unset ztics'
!write ( file_unit, '(a)' ) 'set size ratio 0.5'
!******************************************************************************
write ( file_unit, '(a)' ) 'splot "' // trim ( 'Intensity.txt' ) // &
	& '" using 1:2:3 with lines palette'
write ( file_unit, '(a)' ) 'set view 180, 0'
!write ( file_unit, '(a)' ) 'set key default'
!write ( file_unit, '(a)' ) 'set key left'
write ( file_unit, '(a)' ) 'replot'
!******************************************************************************
!write ( file_unit, '(a)' ) 'splot "' // trim ( data_file_name ) // &
	!& '" wpm3d "contour.txt" w l lc rgb "black"'
!******************************************************************************
if (present(pause)) then
    if (pause<0.0) then
        write ( file_unit, '(a)' ) 'pause -1 "press RETURN to continue"'
        write(*,*)'______________________________________________________________________________'
    else
write(my_pause,'(e9.3)') pause
        write ( file_unit, '(a)' ) 'pause ' // trim(my_pause)
    end if
else
    write ( file_unit, '(a)' ) 'pause 0'
end if
!******************************************************************************
write ( file_unit, '(a)' ) 'q'
close ( unit = file_unit )
!******************************************************************************
call run_gnuplot (command_file_name)
!******************************************************************************
end subroutine twoD
!******************************************************************************
!******************************************************************************
subroutine run_gnuplot(command_file_name)
!***********************************************************************************
implicit none
character (len = 100) command
character (len = *) command_file_name
integer status
integer system
!***********************************************************************************
!  Issue a command to the system that will startup GNUPLOT, using
!  the file we just wrote as input.
!***********************************************************************************
write (command, *) 'gnuplot ' // trim (command_file_name)
status=system(trim(command))
if (status.ne.0) then
    print *,'RUN_GNUPLOT - Fatal error!'
    stop
end if
return
!***********************************************************************************
end subroutine run_gnuplot
!******************************************************************************
!******************************************************************************
!******************************************************************************
subroutine get_unit(iunit)
!***********************************************************************************
implicit none
integer i
integer ios
integer iunit
logical lopen
!***********************************************************************************
iunit=0
do i=1,99
    if (i/= 5 .and. i/=6) then
        inquire (unit=i, opened=lopen, iostat=ios)
        if (ios==0) then
            if (.not.lopen) then
                iunit=i
                return
            end if
            end if
        end if
end do
return
end subroutine get_unit

end module color

!******************************************************************************
!******************************************************************************
     program LDS
     use color
     implicit none
!******************************************************************************
!Declaration of variables
!******************************************************************************
    integer, parameter  :: n1=29, n2=49
     real(kind=8)       :: x, y, AR, V1
     integer            :: i, j, C, DEV
     real               :: d, cmp
     character*10       :: dater
     !character*10		:: dater
!******************************************************************************
     integer      :: IR, IL, VR, VL, NNH, NNV, N, CI(30), disp
     real(kind=8) :: I0, IC(1422,1422), MI(1422,1422), M(1422,1422),MT(1422,1422),Kfact, CR(30)
!******************************************************************************
!generate data for contours
!    do i=1,n1
!        x = i
!    end do
!    do j=1,n2
!        y=j
!        CR=j*i
!    end do
!******************************************************************************
     write(*,*)'                  ______________________________________________________________'
     write(*,*)'                                                                '
     write(*,*)'                            Laboratory Demonstrator System (LDS)'
     write(*,*)'                                                           '
     write(*,*)'                            Geophysics Research Laboratory (GRL)'
     write(*,*)'                               Department of Earth Sciences'
     write(*,*)'                            Eritrea Institute of Technology (EIT)'
     write(*,*)'                                           Eritrea'
     write(*,*)''
     write(*,*)'                                   Copyright (c) 2021 GRL'
     write(*,*)'                                   Under GNU GPL V3 License'
     write(*,*)'                  ______________________________________________________________'
     write(*,*)''
     write(*,*)'ELECTRICAL RESISTIVITY SURVEYING'
     write(*,*)''
!******************************************************************************
     write(*,*) 'Input the date     '
     read(*,*) dater
!Input the required parameters
!******************************************************************************
 10  call InputP(IL, IR, VL, VR, NNH, NNV, d, V1,CI, disp,DEV)
     call Kfactor(IL, IR, VL, VR, Kfact, d)
!******************************************************************************
!Get the matrix of the mesh equations from the programming model
!******************************************************************************
    call Matrix(IL, IR, NNH, NNV, N, M,CI,cmp)
!******************************************************************************
!Get the inverse of the matrix
!******************************************************************************
    call Gaussj(M, MI, N)
!******************************************************************************
!Get the values of the various currents
!******************************************************************************
    call currentV(NNH, NNV, M, V1, IC, I0)
    write(*,*) V1
!******************************************************************************
!Compute the current values of the segments
!******************************************************************************
    call Intensity(IC,NNH,NNV,MT,IL,IR,I0)
!******************************************************************************
!Computing the Apparent resistivity
!******************************************************************************
    call car(VL, VR, IC, I0, Kfact, AR, NNH, NNV, IR, IL, CI, disp, dater, DEV,V1, MT,cmp)
!******************************************************************************
!******************************************************************************
    write(*,*) 'AR     : ', AR
!******************************************************************************
    call twoD(NNV,NNH,CR, pause=-1.0,palette='RGB',terminal='wxt',pm3d='pm3d map impl',persist='no',contour='both',input='')
!******************************************************************************
!Print the output or the result
!******************************************************************************
    write(*,*)'******************************************************************************'
    write(*,*) 'Press 1 for another run or 2-9 to terminate'
    write(*,*)'******************************************************************************'
    read(*,*) C
        if(C.eq.1) goto 10
!******************************************************************************

end program LDS
