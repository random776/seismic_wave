program seismic
    implicit none

    ! 変数宣言
    integer :: ios, i, j, n_out, num
    integer, parameter :: n_angles = 1500, n_points = 1000
    real(8), allocatable :: x(:), y(:), y_old(:), sin_theta(:), cos_theta(:)
    real(8), allocatable :: sin_theta_old(:), theta(:), plus_cos(:)
    real(8), allocatable :: speed_y(:), speed_val(:)
    real(8) :: a, b, t, dt, t_end, pi, y_max, y_min, dy
    character(len=20) :: filename_x, filename_y, filename_t, filename_speed

    ! 配列の確保
    allocate(x(n_angles), y(n_angles), y_old(n_angles), sin_theta(n_angles), &
             cos_theta(n_angles), sin_theta_old(n_angles))
    allocate(theta(n_angles), plus_cos(n_angles))
    allocate(speed_y(n_points), speed_val(n_points))

    ! ---- 初期設定 ----

    ! θの分布（鋭角寄りの分布）
    do i = 1, n_angles
        theta(i) = 90.0 * (1.0_8 - (real(i) / n_angles) ** (1.0_8 / 10.0_8))
    end do

    ! ファイルのオープン（書き出し）
    filename_x = 'seismic_x.dat'
    open(unit=10, iostat=ios, file=trim(filename_x), action='write', form='formatted', status='replace')
    filename_y = 'seismic_y.dat'
    open(unit=11, iostat=ios, file=trim(filename_y), action='write', form='formatted', status='replace')
    filename_t = 'seismic_t.dat'
    open(unit=12, iostat=ios, file=trim(filename_t), action='write', form='formatted', status='replace')
    filename_speed = 'seismic_speed.dat'
    open(unit=13, iostat=ios, file=trim(filename_speed), action='write', form='formatted', status='replace')

    ! 使用する速度プロファイル（1または0）
    num = 1

    ! 初期条件
    x = 0.0_8
    y = 0.0_8
    a = 1.0_8
    b = 2.0_8
    t = 0.0_8
    dt = 0.002_8
    t_end = 10.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8

    ! θからsinとcosを計算
    sin_theta = sin(pi / 180.0_8 * theta)
    cos_theta = cos(pi / 180.0_8 * theta)

    ! ---- 速度プロファイルの出力 ----

    y_max = 10.0_8
    y_min = -10.0_8
    dy = (y_max - y_min) / real(n_points)

    do i = 1, n_points
        speed_y(i) = y_min + dy * i
        speed_val(i) = wave_speed(speed_y(i), num)
        write(13, *) -speed_y(i), speed_val(i)  ! 出力を反転して上向きに
    end do

    ! ---- メインループ：波の追跡 ----

    do while (t < t_end)
        do i = 1, n_angles
            ! 前の状態を保存
            y_old(i) = y(i)
            sin_theta_old(i) = sin_theta(i)

            ! y >= 0 のときのみ更新（媒質内部）
            if (y(i) >= 0.0_8) then
                y(i) = y(i) + cos_theta(i) * wave_speed(y(i), num) * dt
                x(i) = x(i) + sin_theta(i) * wave_speed(y(i), num) * dt
            end if

            ! Snellの法則による屈折の更新
            sin_theta(i) = wave_speed(y(i), num) / wave_speed(y_old(i), num) * sin_theta(i)

            ! 全反射チェック（sin > 1.0）
            if (sin_theta(i) > 1.0_8) then
                sin_theta(i) = sin_theta_old(i)
                plus_cos(i) = -1.0_8  ! 反射したのでcosが反転
            end if

            ! cosの再計算（反射方向も考慮）
            cos_theta(i) = plus_cos(i) * sqrt(1.0_8 - sin_theta(i)**2)
        end do

        ! 時間の更新
        t = t + dt
        n_out = n_out + 1

        ! 一定ステップごとに出力（20ステップごと）
        if (mod(n_out, 20) == 0) then
            write(10, *) x
            write(11, *) -y  ! 描画の都合上、符号を反転
            write(12, *) t
        end if
    end do

    ! ---- 終了処理 ----

    close(10)
    close(11)
    close(12)
    close(13)

    deallocate(x, y, y_old, sin_theta, cos_theta, sin_theta_old)
    deallocate(theta, plus_cos, speed_y, speed_val)

contains

    ! ---- 関数定義：速度プロファイル ----
    real(8) function wave_speed(y, num)
        real(8), intent(in) :: y
        integer, intent(in) :: num

        select case (num)

        ! num = 0: 段階的に速度が増加するプロファイル
        case (0)
            if (0.0_8 <= y .and. y < 4.0_8) then
                wave_speed = y + 1.0_8
            else if (4.0_8 <= y .and. y < 6.0_8) then
                wave_speed = 7.0_8 * y - 23.0_8
            else if (6.0_8 <= y) then
                wave_speed = y + 13.0_8
            else
                wave_speed = 0.0_8
            end if

        ! num = 1: 一度減速してから再加速するプロファイル
        case (1)
            if (0.0_8 <= y .and. y < 4.0_8) then
                wave_speed = y + 1.0_8
            else if (4.0_8 <= y .and. y < 4.375_8) then
                wave_speed = -7.0_8 * y + 33.0_8
            else if (4.375_8 <= y) then
                wave_speed = y - 2.0_8
            else
                wave_speed = 0.0_8
            end if

        ! その他（エラー時）は0を返す
        case default
            wave_speed = 0.0_8
        end select
    end function wave_speed

end program seismic


