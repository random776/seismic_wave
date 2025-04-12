program seismic
    implicit none

    ! ---- 変数宣言 ----
    integer :: ios, i, j, num
    integer, parameter :: n_angles = 500  ! 発射する角度の本数
    real(8), allocatable :: x(:), y(:), y_old(:), sin_theta(:), cos_theta(:)
    real(8), allocatable :: sin_theta_old(:), theta(:), plus_cos(:)
    real(8) :: a, b, t, dt, t_end, pi
    character(len=20) :: filename, filename_curve

    ! ---- 配列の確保 ----
    allocate(x(n_angles), y(n_angles), y_old(n_angles), sin_theta(n_angles), &
             cos_theta(n_angles), sin_theta_old(n_angles))
    allocate(theta(n_angles), plus_cos(n_angles))

    ! ---- 初期角度の設定（鋭角に集中させる）----
    do i = 1, n_angles
        theta(i) = 90.0 * (1.0_8 - (real(i) / n_angles) ** (1.0_8 / 10.0_8))
    end do

    ! ---- 出力ファイルの準備 ----
    filename = 'seismic.dat'  ! 時系列の全体軌跡データ
    open(unit=10, iostat=ios, file=trim(filename), action='write', form='formatted', status='replace')
    filename_curve = 'seismic_curve.dat'  ! 初めて地表に到達した点の記録
    open(unit=11, iostat=ios, file=trim(filename_curve), action='write', form='formatted', status='replace')

    ! ---- シミュレーションパラメータ ----
    num = 1  ! 使用する速度プロファイル（1:一度減速して再加速）
    x = 0.0_8
    y = 0.0_8
    a = 1.0_8  ! 未使用（予備変数？）
    b = 2.0_8  ! 未使用
    t = 0.0_8
    dt = 0.002_8
    t_end = 30.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8  ! 最初はすべて上向きに

    ! ---- sinθ, cosθの初期化 ----
    sin_theta = sin(pi / 180.0_8 * theta)
    cos_theta = cos(pi / 180.0_8 * theta)

    ! ---- メインループ：波線追跡 ----
    do while (t < t_end)
        do i = 1, n_angles
            ! 1ステップ前の情報を保存
            y_old(i) = y(i)
            sin_theta_old(i) = sin_theta(i)

            ! 地表下でのみ位置を更新
            if (y(i) < 0.0_8) then
                y(i) = y(i)
                x(i) = x(i)
            else 
                ! y方向（深さ）とx方向の進行
                y(i) = y(i) + cos_theta(i) * wave_speed(y(i), num) * dt
                x(i) = x(i) + sin_theta(i) * wave_speed(y(i), num) * dt
            end if

            ! 地表（y=0）通過時に位置を記録（上向きに反転）
            if (y(i) * y_old(i) < 0.0_8) then
                write(11, '(F10.4, 2(F15.6, F15.6))') t, x(i), -y(i)
            end if

            ! Snellの法則に基づく角度の変化（速度差を考慮）
            sin_theta(i) = wave_speed(y(i), num) / wave_speed(y_old(i), num) * sin_theta(i)

            ! 全反射が発生した場合
            if (sin_theta(i) > 1.0_8) then
                sin_theta(i) = sin_theta_old(i)  ! 直前の値に戻す
                plus_cos(i) = -1.0_8             ! cos方向を反転（反射）
            end if

            ! 新しいcosθの計算（反射も考慮）
            cos_theta(i) = plus_cos(i) * sqrt(1.0_8 - sin_theta(i)**2)
        end do

        ! 時間更新
        t = t + dt

        ! すべての角度の位置データを1行に出力（行：時刻、列：x, -y, x, -y, ...）
        write(10, '(F10.4, 5(F15.6, F15.6))') t, (x(i), -y(i), i=1,n_angles)
    end do

    ! ---- 後処理 ----
    close(10)
    close(11)

    deallocate(x, y, y_old, sin_theta, cos_theta, sin_theta_old)
    deallocate(theta, plus_cos)

contains

    ! ---- 速度プロファイルを返す関数 ----
    real(8) function wave_speed(y, num)
        real(8), intent(in) :: y
        integer, intent(in) :: num

        select case (num)
        case (0)  ! 単純に深さとともに速度増加
            if (0.0_8 <= y .and. y < 4.0_8) then
                wave_speed = y + 1.0_8
            else if (4.0_8 <= y .and. y < 6.0_8) then
                wave_speed = 7.0_8 * y - 23.0_8
            else if (6.0_8 <= y) then
                wave_speed = y + 13.0_8
            else
                wave_speed = 0.0_8
            end if

        case (1)  ! 一度減速してから再加速
            if (0.0_8 <= y .and. y < 4.0_8) then
                wave_speed = y + 1.0_8
            else if (4.0_8 <= y .and. y < 4.375_8) then
                wave_speed = -7.0_8 * y + 33.0_8
            else if (4.375_8 <= y) then
                wave_speed = y - 2.0_8
            else
                wave_speed = 0.0_8
            end if

        case default  ! 不正な指定
            wave_speed = 0.0_8
        end select
    end function wave_speed

end program seismic

