! high-order simulator for aerodynamics(hosta)

program hosta
    implicit none
    character(len=256) :: msg
    integer(kind=8) :: tm(4), proc_rate, count_max

    call env_initialize

    call system_clock(tm(1), proc_rate, count_max)
    call preset
    call system_clock(tm(2), proc_rate, count_max)
    write(msg, *) "++ preset time cost = ", real(tm(2) - tm(1))/real(proc_rate)
    call msg_seq_and_master(msg)

    call solve
    call system_clock(tm(3), proc_rate, count_max)
    write(msg, *) "++ solver time cost = ", real(tm(3) - tm(2))/real(proc_rate)
    call msg_seq_and_master(msg)

    write(msg, *) "++ total  time cost = ", real(tm(3) - tm(1))/real(proc_rate)
    call msg_seq_and_master(msg)

    call env_finalize

end program hosta
