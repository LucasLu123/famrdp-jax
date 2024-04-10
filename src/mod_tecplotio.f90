
module mod_tecplotio
    implicit none
    private
    public :: tecio_open,tecio_ini,tecio_zone
    public :: tecio_eohmark,tecio_data,tecio_close

    integer,parameter :: int32=4,float32=4
    character(len=1),parameter :: nullchar=char(0)
    character(len=8),parameter :: cmagicnum='#!TDV102'

    integer(int32),parameter :: irworder=1

    real(float32),parameter :: fzonemark=299.0
    real(float32),parameter :: feohmark=357.0

    integer(int32),pointer :: ititle(:)
    integer(int32) :: inumvar
    integer(int32),pointer :: ivarname(:)

    integer(int32),pointer :: izonename(:)
    integer(int32) :: izonecolor
    integer(int32) :: izonetype
    integer(int32) :: idatapack
    integer(int32) :: ivarloc
    integer(int32) :: ifaceconnect
    integer(int32) :: imax,jmax,kmax
    integer(int32) :: iauxiliar

    integer(int32),pointer :: idataform(:)
    integer(int32) :: ivarshare
    integer(int32) :: ishareconnect

    contains

    subroutine tecio_open(fileid,filename)
       implicit none
       integer :: fileid,ierr
       character(len=*) :: filename

       open(fileid,file=filename,form='unformatted',iostat=ierr)

    end subroutine tecio_open

    subroutine tecio_ini(fileid,ctitle,nvar,cvarname)
       implicit none
       character(len=*) :: ctitle
       integer :: fileid,nvar
       character(len=*) :: cvarname
       character(len=len(cvarname)) :: cname
       integer :: i,nchar,nsep,ierr

       write(fileid,iostat=ierr)cmagicnum
       write(fileid,iostat=ierr)irworder

       nchar = len_trim(ctitle)
       allocate(ititle(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,ctitle(1:nchar)//nullchar,ititle(:))
       write(fileid,iostat=ierr)ititle(:)
       deallocate(ititle,stat=ierr)

       inumvar = nvar
       write(fileid,iostat=ierr)inumvar

       nchar = len_trim(cvarname)
       cname = cvarname(1:nchar)
       nsep  = 0
       do i=1,nchar
          if (cname(i:i) == ',') then
             cname(i:i) = nullchar
             nsep = nsep + 1
          end if
       end do
       if (nsep /= inumvar-1) then
          write(*,*) "tecio varname is invalid"
          stop
       end if
       allocate(ivarname(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,cname(1:nchar)//nullchar,ivarname(:))
       write(fileid,iostat=ierr)ivarname(:)
       deallocate(ivarname,stat=ierr)

    end subroutine tecio_ini

    subroutine tecio_zone(fileid,czonename,ndatapack,ni,nj,nk)
       implicit none
       character(len=*) :: czonename
       integer :: fileid,ndatapack,ni,nj,nk
       integer :: nchar,ierr

       write(fileid,iostat=ierr)fzonemark

       nchar = len_trim(czonename)
       allocate(izonename(nchar+1),stat=ierr)
       call tecio_ichar(nchar+1,czonename(1:nchar)//nullchar,izonename(:))
       write(fileid,iostat=ierr)izonename(:)
       deallocate(izonename,stat=ierr)

       izonecolor = -1          !set to -1 if you want tecplot to determine
       write(fileid,iostat=ierr)izonecolor

       izonetype = 0            !0=ORDERED,1=FELINESEG,2=FETRIANGLE,
                                !3=FEQUADRILATERAL,4=FETETRAHEDRON,5=FEBRICK
       write(fileid,iostat=ierr)izonetype

       idatapack = ndatapack    !0=Block, 1=Point
       write(fileid,iostat=ierr)idatapack

       ivarloc = 0              !0 = Don't specify, all data is
                                !located at the nodes. 1 = Specify
       write(fileid,iostat=ierr)ivarloc

       ifaceconnect = 0         !Number of user defined face neighbor connections (value >= 0)
       write(fileid,iostat=ierr)ifaceconnect

       imax = ni
       jmax = nj
       kmax = nk
       write(fileid,iostat=ierr)imax,jmax,kmax

       iauxiliar = 0            !1=Auxiliary name/value pair to follow
                                !0=No more Auxiliar name/value pairs.
       write(fileid,iostat=ierr)iauxiliar

    end subroutine tecio_zone

    subroutine tecio_eohmark(fileid)
       implicit none
       integer :: fileid,ierr

       write(fileid,iostat=ierr)feohmark

    end subroutine tecio_eohmark

    subroutine tecio_data(fileid,ndataform)
       implicit none
       integer :: ndataform
       integer :: fileid,ierr

       write(fileid,iostat=ierr)fzonemark

       allocate(idataform(inumvar),stat=ierr)
       idataform(:) = ndataform !1=Float,2=Double,3=LongInt,4=ShortInt,
                                !5=Byte,6=Bit
       write(fileid,iostat=ierr)idataform(:)
       deallocate(idataform,stat=ierr)

       ivarshare = 0            !0=no,1=yes
       write(fileid,iostat=ierr)ivarshare

       ishareconnect = -1       !Zone number to share connectivity list with (-1 = no sharing).
       write(fileid,iostat=ierr)ishareconnect

    end subroutine tecio_data

    subroutine tecio_close(fileid)
       implicit none
       integer :: fileid,ierr

       close(fileid,iostat=ierr)

    end subroutine tecio_close

    subroutine tecio_ichar(n,ch,ich)
       implicit none
       integer :: n
       character(len=n) :: ch
       integer(int32) :: ich(n)
       integer :: i

       do i=1,n
          ich(i) = ichar(ch(i:i))
       end do

    end subroutine tecio_ichar

end module mod_tecplotio