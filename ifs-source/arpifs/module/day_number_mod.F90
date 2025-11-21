MODULE DAY_NUMBER_MOD
IMPLICIT NONE
PRIVATE
INTEGER, PARAMETER :: N_MONTHS = 12
INTEGER, PARAMETER :: DAYS_PER_MONTH_IN_NONLEAP(N_MONTHS) = &
(/ 31, 28, 31, 30, 31, 30, &
31, 31, 30, 31, 30, 31 /)

PUBLIC :: NUMBER_OF_DAY
CONTAINS

FUNCTION Is_Leap_Year(Year) RESULT(Its_a_Leap_Year)
INTEGER, INTENT(IN) :: Year
LOGICAL :: Its_a_Leap_Year
Its_a_Leap_Year = .FALSE.
IF ( ( MOD( Year, 4 ) == 0 .AND. MOD( Year, 100 ) /= 0 ) .OR. &
MOD( Year, 400 ) == 0 ) Its_a_Leap_Year = .TRUE.
END FUNCTION Is_Leap_Year

FUNCTION Date2DoY( Day_of_Month, & ! Input
Month , & ! Input
Year ) & ! Input
RESULT(Day_Of_Year)
INTEGER, INTENT(IN) :: Day_of_Month
INTEGER, INTENT(IN) :: Month
INTEGER, INTENT(IN) :: Year
INTEGER :: Day_of_Year
INTEGER :: Days_per_Month(N_MONTHS)

Day_Of_Year = -1
! Compute days per month
Days_per_Month = DAYS_PER_MONTH_IN_NONLEAP
IF ( Is_Leap_Year(Year) ) Days_per_Month(2) = 29
! Error checking
IF ( Year < 1 ) RETURN
IF ( Month < 1 .OR. Month > N_MONTHS ) RETURN
IF ( Day_of_Month > Days_per_Month(Month) ) RETURN
! Compute day of year
Day_of_Year = SUM(Days_per_Month(1:Month-1)) + Day_of_Month
END FUNCTION Date2DoY



SUBROUTINE NUMBER_OF_DAY(DD,MM,YYYY,DOY)

INTEGER, INTENT(IN)  :: DD   ! day
INTEGER, INTENT(IN)  :: MM   !month
INTEGER, INTENT(IN)  :: YYYY !year
INTEGER, INTENT(OUT) :: DOY  !number of the day in a certain year (1 to 365 or 366, this is not julian day)

DOY = Date2DoY(DD,MM,YYYY)
RETURN
END SUBROUTINE NUMBER_OF_DAY
 
END MODULE DAY_NUMBER_MOD
