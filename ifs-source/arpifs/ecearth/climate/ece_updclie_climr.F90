subroutine ece_updclie_climr(YDGEOMETRY, YDSURF)

   use PARKIND1, only: JPRB, JPRD, JPIM, JPIB
   use GEOMETRY_MOD, only: GEOMETRY
   use SURFACE_FIELDS_MIX, only: TSURF
   use YOMCT0, only: CNMEXP
   use YOMRIP, only: YRRIP
   use YOMMCC, only: YRMCC
   use YOEPHY, only: YREPHY
   use YOMLUN, only: NULOUT, NULERR
   use YOMMP0, only: MYPROC
   use MPL_MODULE, only: MPL_BROADCAST

   implicit none

   type(GEOMETRY), intent(in) :: YDGEOMETRY
   type(TSURF), intent(inout) :: YDSURF

   integer, parameter :: max_idates = 12*55  ! max no. of dates in the CLIMR file
   integer, parameter :: max_fields = 10  ! max no. of fields (vars) in the file

   type :: t_date
      integer(kind=JPIM) :: ymd  ! YYYYMMDD
      integer(kind=JPIM) :: y, m, d
      integer(kind=JPIM) :: jdn  ! Julian day number (integer part of the Julian date)
   end type

   logical :: first_run = .true.
   integer(kind=JPIM) :: prev = 1, next = 2  ! array indices for previous and next idate
   integer(kind=JPIM) :: log_day = 0  ! used to limit log output to once per day

   integer(kind=JPIM) :: num_fields = 0
   integer(kind=JPIM), save :: field_codes(max_fields)

   integer(kind=JPIM) :: num_idates = 0
   type(t_date), allocatable, save :: idate_toc(:)
   integer(kind=JPIM), save :: first, last
   type(t_date), save :: idate(2)  ! Used to store the prev/next idate

   type(t_date) :: now
   real(KIND=JPRB) :: z_p
   integer(kind=JPIM) :: log_now, broadcast_buf(4)

#include "abor1.intfb.h"
#include "updcal.intfb.h"

   now = current_date()

   if (first_run) then
      call setup_toc()
      if (MYPROC == 1) then
         idate(prev) = prev_idate(now)
         idate(next) = next_idate(now)
         call MPL_BROADCAST( &
            (/idate(1)%ymd, idate(2)%ymd, idate_toc(first)%ymd, idate_toc(last)%ymd/), &
            KTAG=1234321, KROOT=1, CDSTRING='ECE_UPDCLIE_CLIMR')
      else
         call MPL_BROADCAST(broadcast_buf, KTAG=1234321, KROOT=1, CDSTRING='ECE_UPDCLIE_CLIMR')
         idate(1) = date_from_ymd(broadcast_buf(1))
         idate(2) = date_from_ymd(broadcast_buf(2))
         idate_toc(first) = date_from_ymd(broadcast_buf(3))
         idate_toc(last) = date_from_ymd(broadcast_buf(4))
      end if
      call read_fields(idate(prev)%ymd, prev)
      call read_fields(idate(next)%ymd, next)
      first_run = .false.
   end if

   log_now = now%ymd  ! Save current date for logging before shifting
   now = shift_date(now, idate_toc(first), idate_toc(last))

   if (now%ymd == idate(next)%ymd .or. now%ymd < idate(prev)%ymd) then  ! Time to read new fields
      ! Swap values for prev, next: 1 <--> 2
      prev = 3 - prev
      next = 3 - next
      if (MYPROC == 1) then
         idate(prev) = prev_idate(now)  ! prev idate may have changed in case we passed idate_toc(first/last)
         idate(next) = next_idate(now)
         call MPL_BROADCAST(idate%ymd, KTAG=1234322, KROOT=1, CDSTRING='ECE_UPDCLIE_CLIMR')
      else
         call MPL_BROADCAST(broadcast_buf(1:2), KTAG=1234322, KROOT=1, CDSTRING='ECE_UPDCLIE_CLIMR')
         idate = (/date_from_ymd(broadcast_buf(1)), date_from_ymd(broadcast_buf(2))/)
      end if
      call read_fields(idate(next)%ymd, next)
   end if

   ! Compute interpolation weight corresponding to previous idate
   z_p = REAL(idate(next)%jdn - now%jdn, JPRB)/(idate(next)%jdn - idate(prev)%jdn)

   if (MYPROC == 1 .and. now%ymd /= log_day) then
      write (NULOUT, '(a,"[",i8," <",i8," (",i8,")> ",i8,")",2f8.4)') &
         'ece_updclie_climr: Interpolation date (shifted), interval and weights: ', &
         idate(prev)%ymd, log_now, now%ymd, idate(next)%ymd, &
         z_p, 1.0_JPRB - z_p
      log_day = now%ymd
   end if

   call interpolate_fields(prev, z_p)

contains

   pure function jdn(yy, mm, dd)
      ! Julian day number (integer part of Julian date)
      ! https://en.wikipedia.org/wiki/Julian_day#Converting_Gregorian_calendar_date_to_Julian_Day_Number
      integer :: jdn
      integer, intent(in) :: yy, mm, dd
      jdn = (1461*(yy + 4800 + (mm - 14)/12))/4 &
            + (367*(mm - 2 - 12*((mm - 14)/12)))/12 &
            - (3*((yy + 4900 + (mm - 14)/12)/100))/4 &
            + dd - 32075
   end function jdn

   pure function date_from_ymd(ymd) result(date)
      type(t_date) :: date
      integer, intent(in) :: ymd
      date%ymd = ymd
      date%d = MOD(ymd, 100)
      date%m = MOD(ymd/100, 100)
      date%y = ymd/10000
      date%jdn = jdn(date%y, date%m, date%d)
   end function date_from_ymd

   pure function date_from_y_m_d(yy, mm, dd)
      type(t_date) :: date_from_y_m_d
      integer, intent(in) :: yy, mm, dd
      date_from_y_m_d = date_from_ymd(10000*yy + 100*mm + dd)
   end function date_from_y_m_d

   function current_date() result(date)
      use YOMRIP0, only: NINDAT
      use YOMRIP, only: YRRIP
      type(t_date) :: date, start
      integer(kind=JPIM) :: dd, mm, yy, nmm(12)
      start = date_from_ymd(NINDAT)
      call UPDCAL(start%d, start%m, start%y, YRRIP%NSTADD, dd, mm, yy, nmm, -1)
      date = date_from_y_m_d(yy, mm, dd)
   end function current_date

   pure function shift_date(date, date_a, date_b)
      ! Returns a shifted date with same day and months as the date argument,
      ! but within the date interval given by [date_a, date_b) or [date_b, date_a)
      type(t_date) :: shift_date
      type(t_date), intent(in) :: date, date_a, date_b
      type(t_date) :: d_0, d_1
      integer(kind=JPIM) :: yy

      if (date_a%jdn <= date_b%jdn) then
         d_0 = date_a
         d_1 = date_b
      else
         d_0 = date_b
         d_1 = date_a
      end if
      if (d_0%jdn <= date%jdn .and. date%jdn < d_1%jdn) then
         shift_date = date
         return
      end if
      do yy = d_0%y, d_1%y
         shift_date = date_from_y_m_d(yy, date%m, date%d)
         if (d_0%jdn <= shift_date%jdn .and. shift_date%jdn < d_1%jdn) return
      end do
      shift_date = date_from_ymd(0)  ! Panic reaction, could not shift date
   end function shift_date

   pure function prev_idate(date)
      type(t_date) :: prev_idate
      type(t_date), intent(in) :: date
      type(t_date) :: d
      integer :: i

      d = shift_date(date, idate_toc(first), idate_toc(last))
      i = MINLOC(ABS(idate_toc%jdn - d%jdn), dim=1) - 1  ! -1 to account for index 0 in idate_toc
      if (d%jdn < idate_toc(i)%jdn) then
         prev_idate = idate_toc(i - 1)
      else
         prev_idate = idate_toc(i)
      end if
   end function prev_idate

   pure function next_idate(date)
      type(t_date) :: next_idate
      type(t_date), intent(in) :: date
      type(t_date) :: d
      integer :: i

      d = shift_date(date, idate_toc(first), idate_toc(last))
      i = MINLOC(ABS(idate_toc%jdn - d%jdn), dim=1) - 1  ! -1 to account for index 0 in idate_toc
      if (d%jdn < idate_toc(i)%jdn) then
         next_idate = idate_toc(i)
      else
         next_idate = idate_toc(i + 1)
      end if
   end function next_idate

   subroutine setup_toc
      ! Setup idate_toc and num_fields, field_codes
      use GRIB_API_INTERFACE, only: &
         IGRIB_OPEN_FILE, &
         IGRIB_NEW_FROM_FILE, &
         IGRIB_GET_VALUE, &
         IGRIB_RELEASE, &
         IGRIB_CLOSE_FILE, &
         JPGRIB_END_OF_FILE

      type(t_date) :: idates_buf(max_idates)
      character(len=13) :: fname
      integer :: unit, ghandle, rc
      integer :: idt, ifl, i
      integer :: date, gribcode
      type(t_date) :: refdt

      if (MYPROC == 1) then
         fname = 'ICMCL'//CNMEXP(1:4)//'INIT'
         call IGRIB_OPEN_FILE(unit, fname, 'r')
         call IGRIB_NEW_FROM_FILE(unit, ghandle, rc)
         read_dates: do idt = 1, max_idates + 1
            if (idt > max_idates) then
               call ABOR1('ece_updclie_climr (setup_toc): ' &
                          //'Too many dates in ICMCL*INIT. Increase max_idates!')
            end if
            call IGRIB_GET_VALUE(ghandle, 'dataDate', date)
            idates_buf(idt) = date_from_ymd(date)
            read_fields: do ifl = 1, max_fields
               call IGRIB_GET_VALUE(ghandle, 'dataDate', date)
               if (idates_buf(idt)%ymd /= date) exit read_fields
               call IGRIB_GET_VALUE(ghandle, 'paramId', gribcode)
               if (idt == 1) then
                  field_codes(ifl) = gribcode
               else
                  if (field_codes(ifl) /= gribcode) then
                     call ABOR1('ece_updclie_climr (setup_toc): ' &
                                //'Inconsistent order of fields in ICMCL*INIT')
                  end if
               end if
               call IGRIB_RELEASE(ghandle)
               call IGRIB_NEW_FROM_FILE(unit, ghandle, rc)
               if (rc == JPGRIB_END_OF_FILE) exit read_dates
               if (ifl == max_fields) then
                  call ABOR1('ece_updclie_climr (setup_toc): ' &
                             //'Too many fields in ICMCL*INIT. Increase max_fields!')
               end if
            end do read_fields
         end do read_dates
         call IGRIB_CLOSE_FILE(unit)

         num_fields = ifl

         num_idates = idt
         first = 0
         last = num_idates + 1
         allocate (idate_toc(first:last))
         idate_toc(1:num_idates) = idates_buf(1:num_idates)

         ! Set up idate_toc(first)
         refdt = date_from_y_m_d(idate_toc(1)%y + 1, idate_toc(1)%m, idate_toc(1)%d)
         i = 1
         do while (i < num_idates .and. idate_toc(i + 1)%jdn < refdt%jdn)
            i = i + 1
         end do
         idate_toc(first) = date_from_y_m_d(idate_toc(i)%y - 1, idate_toc(i)%m, idate_toc(i)%d)

         ! Set up idate_toc(last)
         refdt = date_from_y_m_d(idate_toc(num_idates)%y - 1, idate_toc(num_idates)%m, idate_toc(num_idates)%d)
         i = num_idates
         do while (i > 1 .and. idate_toc(i - 1)%jdn > refdt%jdn)
            i = i - 1
         end do
         idate_toc(last) = date_from_y_m_d(idate_toc(i)%y + 1, idate_toc(i)%m, idate_toc(i)%d)

         ! Verbose output: idate_toc and field_codes
         write (NULOUT, '(2a,i4,2a)') 'ece_updclie_climr (setup_toc): ', &
            'Found ', num_idates, ' dates in CLIMR file ', fname
         do idt = 1, num_idates
            write (NULOUT, '(a,i4,i10)') &
               'ece_updclie_climr (setup_toc): idate', &
               idt, idate_toc(idt)%ymd
         end do
         write (NULOUT, '(2a,2i10)') 'ece_updclie_climr (setup_toc): ', &
            'Added extra dates before/after idate_toc from CLIMR file: ', &
            idate_toc(first)%ymd, idate_toc(last)%ymd
         write (NULOUT, '(2a,i3,2a)') 'ece_updclie_climr (setup_toc): ', &
            'Found ', num_fields, ' fields in CLIMR file ', fname
         do ifl = 1, num_fields
            write (NULOUT, '(A,i2,i5)') &
               'ece_updclie_climr (setup_toc): field_codes', &
               ifl, field_codes(ifl)
         end do
      else
         first = 1
         last = 2
         allocate (idate_toc(first:last))
      end if

      ! update num_fields on all procs
      call MPL_BROADCAST(num_fields, KTAG=1234323, KROOT=1, CDSTRING='ECE_UPDCLIE_CLIMR')
   end subroutine setup_toc

   subroutine read_fields(date, tlindex)
      use YOMMCC, only: YRMCC
      use DISGRID_MOD, only: DISGRID_SEND, DISGRID_RECV
      use GRIB_API_INTERFACE, only: &
         IGRIB_OPEN_FILE, &
         IGRIB_NEW_FROM_FILE, &
         IGRIB_GET_VALUE, &
         IGRIB_RELEASE, &
         IGRIB_CLOSE_FILE, &
         JPGRIB_END_OF_FILE, &
         JPGRIB_SUCCESS

      integer(kind=JPIM), intent(in) :: date
      integer, intent(in) :: tlindex  ! time level index (1 or 2)

      logical :: first_run = .true.
      integer, save :: unit, ghandle
      character(len=13), save :: fname

      integer(kind=JPIM) :: read_date, d, pid, rc, i
      real(kind=JPRB) :: glo_fieldbuf(YDGEOMETRY%YRGEM%NGPTOTG)
      real(kind=JPRB) :: loc_fieldbuf(YDGEOMETRY%YRGEM%NGPTOT)

      if (MYPROC == 1) then
         ! Adapt date for first/last entry in idate_toc
         if (date == idate_toc(first)%ymd) then
            read_date = idate_toc(first)%ymd + 10000
         elseif (date == idate_toc(last)%ymd) then
            read_date = idate_toc(last)%ymd - 10000
         else
            read_date = date
         end if

         write (NULOUT, '(2a,i8)') 'ece_updclie_climr (read_fields): ', &
            'Reading CLIMR fields for date ', read_date

         if (first_run) then  ! Open grib file on first call
            fname = 'ICMCL'//CNMEXP(1:4)//'INIT'
            write (NULOUT, '(3a)') 'ece_updclie_climr (read_fields): ', &
               'Opening CLIMR file ', fname
            call IGRIB_OPEN_FILE(unit, fname, 'r')
            call IGRIB_NEW_FROM_FILE(unit, ghandle, rc)
            if (rc /= JPGRIB_SUCCESS) then
               call ABOR1('ece_updclie_climr (read_fields): Could not open CLIMR file for reading')
            end if
            first_run = .false.
         end if
         call IGRIB_GET_VALUE(ghandle, 'dataDate', d)
         if (d > read_date) then  ! Rewind file
            write (NULOUT, '(3a)') 'ece_updclie_climr (read_fields): ', &
               'Rewind CLIMR file ', fname
            call IGRIB_CLOSE_FILE(unit)
            call IGRIB_OPEN_FILE(unit, fname, 'r')
            call IGRIB_NEW_FROM_FILE(unit, ghandle)
            call IGRIB_GET_VALUE(ghandle, 'dataDate', d, rc)
         end if
         if (d < read_date) then  ! Forward
            do while (d < read_date)
               call IGRIB_RELEASE(ghandle)
               call IGRIB_NEW_FROM_FILE(unit, ghandle, rc)
               if (rc == JPGRIB_END_OF_FILE) then
                  call ABOR1( &
                     'ece_updclie_climr (read_fields): ' &
                     //'Could not find requested date in file')
               end if
               call IGRIB_GET_VALUE(ghandle, 'dataDate', d)
            end do
         end if
         if (d /= read_date) then  ! Date not found
            write (NULOUT, '(a,i8,a)') &
               'ece_updclie_climr (read_fields): Trying to locate date ', &
               read_date, ' in CLIMR file but could not find it'
            call ABOR1('ece_updclie_climr (read_fields): Could not find date in CLIMR file')
         end if
      end if

      ! Found date, now read and distribute
      do i = 1, num_fields
         if (MYPROC == 1) then
            call IGRIB_GET_VALUE(ghandle, 'dataDate', d)
            call IGRIB_GET_VALUE(ghandle, 'paramId', pid)
            call IGRIB_GET_VALUE(ghandle, 'values', glo_fieldbuf)

            if (d /= read_date) then
               write (NULERR, '(2a,i4,2i10)') 'ece_updclie_climr (read_fields): ', &
                  'Date mismatch: field no., expected date, actual: ', &
                  i, read_date, d
               call ABOR1('ece_updclie_climr (read_fields): Date mismatch')
            else if (pid /= field_codes(i)) then
               write (NULERR, '(2a,i4,i10,2i5)') 'ece_updclie_climr (read_fields): ', &
                  'GRIB code mismatch: field no., date, expected code, actual: ', &
                  i, d, field_codes(i), pid
               call ABOR1('ece_updclie_climr (read_fields): GRIB code mismatch')
            else
               write (NULOUT, '(2a,i4,i10,i5)') 'ece_updclie_climr (read_fields): ', &
                  'Read CLIMR field: no., date, gribcode:', i, d, pid
            end if

            call DISGRID_SEND(YDGEOMETRY, 1, glo_fieldbuf, i, loc_fieldbuf)

            if (i < num_fields) then
               call IGRIB_RELEASE(ghandle)
               call IGRIB_NEW_FROM_FILE(unit, ghandle, rc)
               if (rc == JPGRIB_END_OF_FILE .and. i < num_fields) then
                  call ABOR1('ece_updclie_climr (read_fields): ' &
                             //'End-of-file before all fields could be read')
               else if (rc /= JPGRIB_SUCCESS) then
                  call ABOR1('ece_updclie_climr (read_fields): ' &
                             //'Error while trying to read next field from ICMCL file')
               end if
            end if
         else
            call DISGRID_RECV(YDGEOMETRY%YRGEM, 1, 1, loc_fieldbuf, i)
         end if
         YRMCC%CLIMR(:, tlindex, i) = REAL(loc_fieldbuf(:), JPRD)
      end do
   end subroutine read_fields

   subroutine interpolate_fields(np, zp)
      integer(kind=JPIM), intent(in) :: np  ! Index for previous time level (1 or 2)
      real(kind=JPRB), intent(in) :: zp  ! Interpolation weight for previous time level

      integer(kind=JPIM) :: jstglo, jrof, iend, ibl
      integer(kind=JPIM) :: gc66, gc67, gc15, gc16, gc17, gc18, gc174  ! Indices for fields with respective GRIB codes
      integer(kind=JPIM) :: nn  ! Index for next time level (1 or 2)
      real(kind=JPRB) :: zn  ! Interpolation weight for next time level

      associate ( &
         NGPTOT => YDGEOMETRY%YRGEM%NGPTOT, &
         NPROMA => YDGEOMETRY%YRDIM%NPROMA, &
         LE4ALB => YREPHY%LE4ALB, &
         NCLIGC => YRMCC%NCLIGC, &
         SD_VF => YDSURF%SD_VF, &
         YSD_VF => YDSURF%YSD_VF, &
         CLIMR => YRMCC%CLIMR &
         )

         nn = 3 - np
         zn = 1.0_JPRB - zp

         gc174 = findloc(NCLIGC, 174)
         gc15 = findloc(NCLIGC, 15)
         gc16 = findloc(NCLIGC, 16)
         gc17 = findloc(NCLIGC, 17)
         gc18 = findloc(NCLIGC, 18)
         gc66 = findloc(NCLIGC, 66)
         gc67 = findloc(NCLIGC, 67)

         do jstglo = 1, NGPTOT, NPROMA
            iend = MIN(NPROMA, NGPTOT - jstglo + 1)
            ibl = (jstglo - 1)/NPROMA + 1

            do jrof = 1, iend

               ! Albedo, grib parameter 174
               SD_VF(jrof, YSD_VF%YALBF%MP, ibl) = &
                  zp*CLIMR(jstglo + jrof - 1, np, gc174) &
                  + zn*CLIMR(jstglo + jrof - 1, nn, gc174)

               ! MODIS albedo, grib parameters 15, 16, 17, 18
               if (LE4ALB) then
                  SD_VF(jrof, YSD_VF%YALUVP%MP, ibl) = &
                     zp*CLIMR(jstglo + jrof - 1, np, gc15) &
                     + zn*CLIMR(jstglo + jrof - 1, nn, gc15)
                  SD_VF(jrof, YSD_VF%YALUVD%MP, ibl) = &
                     zp*CLIMR(jstglo + jrof - 1, np, gc16) &
                     + zn*CLIMR(jstglo + jrof - 1, nn, gc16)
                  SD_VF(jrof, YSD_VF%YALNIP%MP, ibl) = &
                     zp*CLIMR(jstglo + jrof - 1, np, gc17) &
                     + zn*CLIMR(jstglo + jrof - 1, nn, gc17)
                  SD_VF(jrof, YSD_VF%YALNID%MP, ibl) = &
                     zp*CLIMR(jstglo + jrof - 1, np, gc18) &
                     + zn*CLIMR(jstglo + jrof - 1, nn, gc18)
               end if

               ! LAI low/high, grib parameters 66, 67
               SD_VF(jrof, YSD_VF%YLAIL%MP, ibl) = &
                  zp*CLIMR(jstglo + jrof - 1, np, gc66) &
                  + zn*CLIMR(jstglo + JROF - 1, nn, gc66)
               SD_VF(jrof, YSD_VF%YLAIH%MP, ibl) = &
                  zp*CLIMR(jstglo + jrof - 1, np, gc67) &
                  + zn*CLIMR(jstglo + jrof - 1, nn, gc67)

            end do
         end do

      end associate
   end subroutine interpolate_fields

   function findloc(array, value)
      ! Partial implementation of Fortran findloc intrinsic (Fortran2008)
      integer(kind=JPIM) :: findloc
      integer(kind=JPIM), intent(in) :: array(:), value
      integer :: i
      do i = 1, UBOUND(array, 1)
         if (array(i) == value) exit
      end do
      if (i <= UBOUND(array, 1)) then
         findloc = i
      else
         findloc = 0
      end if
   end function findloc

end subroutine ece_updclie_climr
