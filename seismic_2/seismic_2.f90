program seismic
    implicit none
    integer :: ios, i, j
    integer, parameter :: n_angles = 500
    real(8), allocatable :: x(:), y(:), y_old(:), sin_theta(:), cos_theta(:)
    real(8), allocatable :: sin_theta_old(:), theta(:), plus_cos(:)
    real(8) :: a, b, t, dt, t_end, pi
    character(len=20) :: filename, filename_curve

    allocate(x(n_angles), y(n_angles), y_old(n_angles), sin_theta(n_angles), &
    cos_theta(n_angles), sin_theta_old(n_angles))
    allocate(theta(n_angles), plus_cos(n_angles))

    do i = 1, n_angles
        theta(i) = 90.0 * (1.0_8 - (real(i) / n_angles) ** (1.0_8 / 10.0_8))
    end do

    ! 結果を1つのファイルにまとめて保存
    filename = 'seismic.dat'
    open(unit=10, iostat=ios, file=trim(filename), action='write', form='formatted', status='replace')
    filename_curve = 'seismic_curve.dat'
    open(unit=11, iostat=ios, file=trim(filename_curve), action='write', form='formatted', status='replace')

    ! 初期条件
    x = 0.0_8
    y = 0.0_8
    a = 1.0_8
    b = 2.0_8
    t = 0.0_8
    dt = 0.002_8
    t_end = 30.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8

    sin_theta = sin(pi / 180.0_8 * theta)
    cos_theta = cos(pi / 180.0_8 * theta)


    do while (t < t_end)
        do i = 1, n_angles
            y_old(i) = y(i)
            sin_theta_old(i) = sin_theta(i)
            if (y(i) < 0.0_8) then
                y(i) = y(i)
                x(i) = x(i)
            else 
                y(i) = y(i) + cos_theta(i) * wave_speed(y(i)) * dt
                x(i) = x(i) + sin_theta(i) * wave_speed(y(i)) * dt
            end if
            if (y(i) * y_old(i) < 0.0_8) then
                write(11, '(F10.4, 2(F15.6, F15.6))') t, x(i), -y(i)
            end if
            sin_theta(i) = wave_speed(y(i)) / wave_speed(y_old(i)) * sin_theta(i)
            if (sin_theta(i) > 1.0_8) then
                sin_theta(i) = sin_theta_old(i)
                plus_cos(i) = - 1.0_8
            end if
            cos_theta(i) = plus_cos(i) * sqrt(1.0_8 - sin_theta(i)**2)
        end do
        t = t + dt

        ! すべての角度のデータを1行にまとめて書き込み
        write(10, '(F10.4, 5(F15.6, F15.6))') t, (x(i), -y(i), i=1,n_angles)
    end do

    ! ファイルを閉じる
    close(10)
    close(11)

    deallocate(x, y, y_old, sin_theta, cos_theta, sin_theta_old)
    deallocate(theta, plus_cos)

contains

real(8) function wave_speed(y)
    real(8), intent(in) :: y
    if (0.0_8 <= y .and. y < 4.0_8) then
        wave_speed = y + 1.0_8
    else if (4.0_8 <= y .and. y < 6.0_8) then
        wave_speed = 7.0_8 * y - 23.0_8
    else if (6.0_8 <= y) then
        wave_speed = y + 13.0_8
    else
        wave_speed = 0.0_8
    end if
end function wave_speed

end program seismic
