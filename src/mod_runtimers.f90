
module mod_runtimers
    use mod_kndconsts, only : kind_int,kind_real
    implicit none
    private
    integer(kind_int) :: count_start,count_end,diff_count
    integer(kind_int) :: proc_rate,count_max
    real(kind_real)   :: time_wall,time_wall_sav
    real(kind_real)   :: time_begin,time_end,diff_time
    real(kind_real)   :: time_cpu,time_cpu_sav
    integer(kind_int) :: time_hour,time_minute,time_second
    character(len=12) :: ctime_cpu,ctime_wall

    public :: proc_rate
    public :: time_begin,time_end
    public :: time_cpu,time_wall
    public :: time_cpu_sav,time_wall_sav
    public :: start_timer,get_run_timer
    public :: ctime_cpu,ctime_wall

    contains

    subroutine start_timer
       implicit none

       call start_cpu_timer
       call start_wall_timer

    end subroutine start_timer

    subroutine get_run_timer
       implicit none

       call get_cpu_timer
       call get_wall_timer

    end subroutine get_run_timer

    subroutine start_wall_timer
       implicit none

       call system_clock(count_start, proc_rate, count_max)

       time_wall = time_wall_sav

    end subroutine start_wall_timer

    subroutine get_wall_timer
       implicit none

       call system_clock(count_end)

       diff_count = count_end - count_start
       if (diff_count < 0) then
          diff_count = count_max + diff_count
       end if
       count_start = count_end

       time_wall = time_wall + diff_count/real(proc_rate)

       call second_to_hour_minute(time_wall)
       call sav_time_to_character(ctime_wall)

    end subroutine get_wall_timer

    subroutine start_cpu_timer
        implicit none

        call cpu_time(time_begin)

    end subroutine start_cpu_timer

    subroutine get_cpu_timer
        implicit none

        call cpu_time(time_end)

        time_cpu = time_cpu_sav

        diff_time = time_end - time_begin

        time_cpu = time_cpu + diff_time

        call second_to_hour_minute(time_cpu)
        call sav_time_to_character(ctime_cpu)

    end subroutine get_cpu_timer

    subroutine second_to_hour_minute(time)
        implicit none
        real(kind_real) :: time

        time_hour   = int(time/3600.0)
        time_minute = int(mod(time,3600.0)/60.0)
        time_second = int(mod(mod(time,3600.0),60.0))

    end subroutine second_to_hour_minute

    subroutine sav_time_to_character(ctime)
        implicit none
        character(len=*) :: ctime

        if (time_minute < 10) then
            if (time_second < 10) then
                 write(ctime,'(i5,":0",i1,":0",i1)')time_hour,time_minute,time_second
            else
                 write(ctime,'(i5,":0",i1,":",i2)')time_hour,time_minute,time_second
            end if
        else
            if (time_second < 10) then
                write(ctime,'(i5,":",i2,":0",i1)')time_hour,time_minute,time_second
            else
                write(ctime,'(i5,":",i2,":",i2)')time_hour,time_minute,time_second
            end if
        end if

    end subroutine sav_time_to_character

end module mod_runtimers
