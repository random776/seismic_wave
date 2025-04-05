program seismic
    implicit none
    integer :: ios
    real(8) :: x, y, y_old, a, b, sin_theta, cos_theta, sin_theta_old, theta, t, dt, t_end, pi, plus_cos

    open(unit=10, iostat=ios, file='seismic_1.dat', action='write', form='formatted', status='replace')
    open(unit=11, iostat=ios, file='seismic_2.dat', action='write', form='formatted', status='replace')
    open(unit=12, iostat=ios, file='seismic_3.dat', action='write', form='formatted', status='replace')

    x = 0.0_8
    y = 0.0_8
    a = 1.0_8
    b = 2.0_8
    theta = 30.0_8
    t = 0.0_8
    dt = 0.01_8
    t_end = 10.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8

    sin_theta = sin(pi / 180.0_8 * theta)
    cos_theta = cos(pi / 180.0_8 * theta)

    write(*, *) 'Initial cos: ', pi / 180.0_8 * theta


    do while (t < t_end .and. y >= 0.0_8)
        y_old = y
        sin_theta_old = sin_theta
        y = y + cos_theta * wave_speed(y, a, b) * dt ! x, yを1ステップ進める
        x = x + sin_theta * wave_speed(y, a, b) * dt 
        sin_theta = wave_speed(y, a, b) / wave_speed(y_old, a, b) * sin_theta ! スネルの法則
        if (sin_theta > 1.0_8) then ! 降下→上昇の切り替え
            sin_theta = sin_theta_old
            plus_cos = - 1.0_8
        end if
        cos_theta = plus_cos * sqrt(1.0_8 - sin_theta**2)
        t = t + dt
        write(10, '(F20.12)') t
        write(11, '(F20.12)') x
        write(12, '(F20.12)') - y
    end do

    close(10)

contains

real(8) function wave_speed(y, a, b) ! 速度プロファイル
        real(8) :: y, a, b

        wave_speed = a * y + b

    end function wave_speed

end program seismic