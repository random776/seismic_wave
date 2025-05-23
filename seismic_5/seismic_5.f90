program seismic
    implicit none
    ! 変数の宣言
    integer :: ios, i, j, n_out
    integer, parameter :: n_angles = 5000, n_points = 100
    real(8), allocatable :: x(:), y(:), x_old(:), y_old(:), sin_theta(:), cos_theta(:)
    real(8), allocatable :: sin_theta_old(:), theta(:), plus_cos(:)
    real(8), allocatable :: speed_x(:), speed_y(:), speed_val(:, :)
    real(8), allocatable :: sin_theta_mod(:), cos_theta_mod(:), x_mod(:), y_mod(:)
    real(8) :: a, b, t, dt, div, t_end, pi, x_max, x_min, dx, y_max, y_min, dy, x_prev, y_prev
    real(8) :: delv_delx, delv_dely
    character(len=20) :: filename_x, filename_y, filename_t, filename_speed, filename_curve

    ! 配列の割り当て
    allocate(x(n_angles), y(n_angles), x_old(n_angles), y_old(n_angles))
    allocate(sin_theta(n_angles), cos_theta(n_angles), sin_theta_old(n_angles))
    allocate(theta(n_angles), plus_cos(n_angles))
    allocate(speed_x(n_points), speed_y(n_points), speed_val(n_points, n_points))
    allocate(sin_theta_mod(n_angles), cos_theta_mod(n_angles), x_mod(n_angles), y_mod(n_angles))

    ! 入射角の設定（浅い角度を多くサンプリング）
    do i = 1, n_angles
        theta(i) = - 90.0 * (real(i) / n_angles) ** (1.0_8 / 15.0_8) !負であることに注意
    end do

    ! 出力ファイルの準備
    filename_x = 'seismic_x.dat'
    open(unit=10, iostat=ios, file=trim(filename_x), action='write', form='formatted', status='replace')
    filename_y = 'seismic_y.dat'
    open(unit=11, iostat=ios, file=trim(filename_y), action='write', form='formatted', status='replace')
    filename_t = 'seismic_t.dat'
    open(unit=12, iostat=ios, file=trim(filename_t), action='write', form='formatted', status='replace')
    filename_speed = 'seismic_speed.dat'
    open(unit=13, iostat=ios, file=trim(filename_speed), action='write', form='formatted', status='replace')
    filename_curve = 'seismic_curve.dat'
    open(unit=14, iostat=ios, file=trim(filename_curve), action='write', form='formatted', status='replace')

    ! 初期設定
    x = 0.02_8
    y = 9.98_8         ! 発信点の初期高さ（地表）
    a = 1.0_8
    b = 2.0_8
    t = 0.0_8
    dt = 0.001_8
    div = 0.01_8
    t_end = 6.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8

    theta = theta * pi / 180.0_8

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    ! write(*, *) theta

    ! 波速の分布（グリッドデータ）を生成
    x_max = 20.0_8
    x_min = -20.0_8
    y_max = 20.0_8
    y_min = -20.0_8
    dx = (x_max-x_min) / real(n_points)
    dy = (y_max-y_min) / real(n_points)
    do i = 1, n_points
        do j = 1, n_points
            speed_x(i) = x_min + dx * i
            speed_y(j) = y_min + dy * j
            speed_val(i, j) = wave_speed(speed_x(i), speed_y(j))
            write(13, *) speed_x(i), speed_y(j), speed_val(i, j)
        end do
    end do

    ! メインの時間発展ループ
    do while (t < t_end)
        do i = 1, n_angles
            ! 伝播中かどうかを判定
            if (sqrt(x(i)**2 + y(i)**2) <= 10.0_8 .and. x(i) >= 0.0_8 .and. y(i) >= 0.0_8) then

                ! 位置の更新
                x_prev = x(i)
                y_prev = y(i)
                
                delv_delx = (wave_speed(x(i) + div, y(i)) - wave_speed(x(i) - div, y(i))) / (2*div)
                delv_dely = (wave_speed(x(i), y(i) + div) - wave_speed(x(i), y(i) - div)) / (2*div)

                theta(i) = theta(i) + dt * (delv_delx * sin_theta(i) - delv_dely * cos_theta(i))
                sin_theta(i) = sin(theta(i))
                cos_theta(i) = cos(theta(i))
                x(i) = x(i) + dt * wave_speed(x(i), y(i)) * cos_theta(i)
                y(i) = y(i) + dt * wave_speed(x(i), y(i)) * sin_theta(i)

                ! 境界到達時に曲線情報を記録
                if (sqrt(x_prev ** 2 + y_prev ** 2) <= 10.0_8 .and. sqrt(x(i) ** 2 + y(i) ** 2) > 10.0_8) then
                    write(14, *) t, sqrt(x(i) ** 2 + y(i) ** 2) * acos(y(i) / sqrt(x(i) ** 2 + y(i) ** 2))
                end if
            end if
        end do

        ! 時間更新
        t = t + dt
        n_out = n_out + 1

        ! 一定ステップごとに出力
        if (mod(n_out, 50) == 0) then
            write(10, *) x
            write(11, *) y
            write(12, *) t
        end if
    end do

    ! ファイルのクローズ
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)

    ! メモリの解放
    deallocate(x, y, x_mod, y_mod, x_old, y_old, sin_theta, cos_theta, sin_theta_mod, cos_theta_mod)
    deallocate(sin_theta_old, speed_x, speed_y, speed_val)

contains

! 半径に応じた波の伝播速度
real(8) function wave_speed(x, y)
    real(8), intent(in) :: x, y
    real(8) :: r
    r = sqrt(x**2 + y**2)
    if (7.45_8 <= r .and. r <= 10.0_8) then
        wave_speed = 11.0_8 - r
    else if (6.5_8 <= r .and. r < 7.45_8) then
        wave_speed = 160.0_8 - 21.0_8 * r
    else if (0.0_8 <= r .and. r < 6.5_8) then
        wave_speed = 30.0_8 - r
    else
        wave_speed = 0.0_8
    end if
end function wave_speed

end program seismic