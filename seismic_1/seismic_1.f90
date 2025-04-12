program seismic
    implicit none

    ! 変数の宣言
    integer :: ios                      ! 入出力ステータス
    real(8) :: x, y, y_old              ! 現在・過去の座標
    real(8) :: a, b                     ! 速度プロファイルのパラメータ
    real(8) :: sin_theta, cos_theta     ! sinθ, cosθ の現在値
    real(8) :: sin_theta_old            ! 前のステップの sinθ
    real(8) :: theta                    ! 初期角度 [deg]
    real(8) :: t, dt, t_end             ! 時間とタイムステップ
    real(8) :: pi                       ! 円周率
    real(8) :: plus_cos                 ! cos の符号制御（下向き or 上向き）

    ! 出力ファイルを開く
    open(unit=10, iostat=ios, file='seismic_1.dat', action='write', form='formatted', status='replace') ! 時間
    open(unit=11, iostat=ios, file='seismic_2.dat', action='write', form='formatted', status='replace') ! x座標
    open(unit=12, iostat=ios, file='seismic_3.dat', action='write', form='formatted', status='replace') ! y座標

    ! 初期条件の設定
    x = 0.0_8
    y = 0.0_8
    a = 1.0_8                 ! y依存の速度傾き
    b = 2.0_8                 ! y=0での初期波速度
    theta = 30.0_8           ! 発射角度 [度]
    t = 0.0_8
    dt = 0.01_8
    t_end = 10.0_8
    pi = 4.0_8 * atan(1.0_8)
    plus_cos = 1.0_8         ! 最初は下向きに進行

    ! 角度をラジアンに変換
    sin_theta = sin(pi / 180.0_8 * theta)
    cos_theta = cos(pi / 180.0_8 * theta)

    write(*, *) 'Initial cos: ', pi / 180.0_8 * theta

    ! 波の進行ループ（地表に戻るまで or 時間終了まで）
    do while (t < t_end .and. y >= 0.0_8)
        y_old = y
        sin_theta_old = sin_theta

        ! y方向に1ステップ進める（cosθに対応）
        y = y + cos_theta * wave_speed(y, a, b) * dt

        ! x方向に1ステップ進める（sinθに対応）
        x = x + sin_theta * wave_speed(y, a, b) * dt

        ! スネルの法則に従って sinθ を更新
        sin_theta = wave_speed(y, a, b) / wave_speed(y_old, a, b) * sin_theta

        ! sinθが1を超えると全反射：符号反転して上向き進行
        if (sin_theta > 1.0_8) then
            sin_theta = sin_theta_old     ! 直前の値に戻す
            plus_cos = -1.0_8             ! cosθの符号を反転（上向きへ）
        end if

        ! cosθ を三角関数の恒等式で更新
        cos_theta = plus_cos * sqrt(1.0_8 - sin_theta**2)

        ! 時間を更新
        t = t + dt

        ! 出力：時刻、x座標、y座標（マイナスにして下が正方向になるように）
        write(10, '(F20.12)') t
        write(11, '(F20.12)') x
        write(12, '(F20.12)') - y
    end do

    ! ファイルを閉じる
    close(10)
    close(11)
    close(12)

contains

    !---------------------------------------------------
    ! 任意のyにおける波の速度：線形に増加するプロファイル
    ! v(y) = a * y + b
    !---------------------------------------------------
    real(8) function wave_speed(y, a, b)
        real(8) :: y, a, b
        wave_speed = a * y + b
    end function wave_speed

end program seismic