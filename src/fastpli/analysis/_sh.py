import numpy as np
from numba import njit


@njit(cache=True)
def _spherical_harmonics_band0(
        order: int, costheta: np.ndarray, sintheta: np.ndarray,
        phi: np.ndarray) -> np.float32:  # pylint: disable=unused-argument
    return 0.2820947917738781434740397257803862929220 + 0 * phi


@njit(cache=True)
def _spherical_harmonics_band2(order: int, costheta: np.ndarray,
                               sintheta: np.ndarray,
                               phi: np.ndarray) -> np.float32:
    if order == -2:
        return 0.5462742152960395352716928529013442013451 * sintheta**2 * np.sin(
            2 * phi)
    elif order == -1:
        return -1.0925484305920790705433857058026884026900 * sintheta * costheta * np.sin(
            phi)
    elif order == 0:
        return 0.3153915652525200060308936902957104933242 * (3 * costheta**2 -
                                                             1)
    elif order == 1:
        return -1.0925484305920790705433857058026884026900 * sintheta * costheta * np.cos(
            phi)
    elif order == 2:
        return 0.5462742152960395352716928529013442013451 * sintheta**2 * np.cos(
            2 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def _spherical_harmonics_band4(order: int, costheta: np.ndarray,
                               sintheta: np.ndarray,
                               phi: np.ndarray) -> np.float32:
    if order == -4:
        return 0.6258357354491761345866405236050888509857 * sintheta**4 * np.sin(
            4 * phi)
    elif order == -3:
        return -1.7701307697799305310368308326244860500080 * sintheta**3 * costheta * np.sin(
            3 * phi)
    elif order == -2:
        return 0.4730873478787800090463405354435657399864 * sintheta**2 * (
            7 * costheta**2 - 1) * np.sin(2 * phi)
    elif order == -1:
        return -0.6690465435572891679521123897119059713981 * sintheta * (
            7 * costheta**3 - 3 * costheta) * np.sin(phi)
    elif order == 0:
        return 0.1057855469152043038027648971676448598457 * (
            35 * costheta**4 - 30 * costheta**2 + 3)
    elif order == 1:
        return -0.6690465435572891679521123897119059713981 * sintheta * (
            7 * costheta**3 - 3 * costheta) * np.cos(phi)
    elif order == 2:
        return 0.4730873478787800090463405354435657399864 * sintheta**2 * (
            7 * costheta**2 - 1) * np.cos(2 * phi)
    elif order == 3:
        return -1.7701307697799305310368308326244860500080 * sintheta**3 * costheta * np.cos(
            3 * phi)
    elif order == 4:
        return 0.6258357354491761345866405236050888509857 * sintheta**4 * np.cos(
            4 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def _spherical_harmonics_band6(order: int, costheta: np.ndarray,
                               sintheta: np.ndarray,
                               phi: np.ndarray) -> np.float32:
    if order == -6:
        return 0.6831841051919143219748071752758860124128 * sintheta**6 * np.sin(
            6 * phi)
    elif order == -5:
        return -2.3666191622317520319877380468747716155400 * sintheta**5 * costheta * np.sin(
            5 * phi)
    elif order == -4:
        return 0.5045649007287241592544788598006009966061 * sintheta**4 * (
            11 * costheta**2 - 1) * np.sin(4 * phi)
    elif order == -3:
        return -0.9212052595149234991415657962237370695259 * sintheta**3 * (
            11 * costheta**3 - 3 * costheta) * np.sin(3 * phi)
    elif order == -2:
        return 0.4606026297574617495707828981118685347630 * sintheta**2 * (
            33 * costheta**4 - 18 * costheta**2 + 1) * np.sin(2 * phi)
    elif order == -1:
        return -0.5826213625187313888350785769893992527669 * sintheta * (
            33 * costheta**5 - 30 * costheta**3 + 5 * costheta) * np.sin(phi)
    elif order == 0:
        return 0.0635692022676284259328270310605563631508 * (
            231 * costheta**6 - 315 * costheta**4 + 105 * costheta**2 - 5)
    elif order == 1:
        return -0.5826213625187313888350785769893992527669 * sintheta * (
            33 * costheta**5 - 30 * costheta**3 + 5 * costheta) * np.cos(phi)
    elif order == 2:
        return 0.4606026297574617495707828981118685347630 * sintheta**2 * (
            33 * costheta**4 - 18 * costheta**2 + 1) * np.cos(2 * phi)
    elif order == 3:
        return -0.9212052595149234991415657962237370695259 * sintheta**3 * (
            11 * costheta**3 - 3 * costheta) * np.cos(3 * phi)
    elif order == 4:
        return 0.5045649007287241592544788598006009966061 * sintheta**4 * (
            11 * costheta**2 - 1) * np.cos(4 * phi)
    elif order == 5:
        return -2.3666191622317520319877380468747716155400 * sintheta**5 * costheta * np.cos(
            5 * phi)
    elif order == 6:
        return 0.6831841051919143219748071752758860124128 * sintheta**6 * np.cos(
            6 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def _spherical_harmonics_band8(order: int, costheta: np.ndarray,
                               sintheta: np.ndarray,
                               phi: np.ndarray) -> np.float32:
    if order == -8:
        return 0.7289266601748298688704348855571033335207 * sintheta**8 * np.sin(
            8 * phi)
    elif order == -7:
        return -2.9157066406993194754817395422284133340830 * sintheta**7 * costheta * np.sin(
            7 * phi)
    elif order == -6:
        return 0.5323327660595425743111699467876298815112 * sintheta**6 * (
            15 * costheta**2 - 1) * np.sin(6 * phi)
    elif order == -5:
        return -3.4499106220981080175178857353706008701770 * sintheta**5 * (
            5 * costheta**3 - costheta) * np.sin(5 * phi)
    elif order == -4:
        return 0.4784165247593306972508126646608462146527 * sintheta**4 * (
            65 * costheta**4 - 26 * costheta**2 + 1) * np.sin(4 * phi)
    elif order == -3:
        return -1.2352661552955440759618898084529796356280 * sintheta**3 * (
            39 * costheta**5 - 26 * costheta**3 + 3 * costheta) * np.sin(
                3 * phi)
    elif order == -2:
        return 0.4561522584349094705715788531131085781609 * sintheta**2 * (
            143 * costheta**6 - 143 * costheta**4 + 33 * costheta**2 -
            1) * np.sin(2 * phi)
    elif order == -1:
        return -0.1090412458987799555260481189135496447649 * sintheta * (
            715 * costheta**7 - 1001 * costheta**5 + 385 * costheta**3 -
            35 * costheta) * np.sin(phi)
    elif order == 0:
        return 0.0090867704915649962938373432427958037304 * (
            6435 * costheta**8 - 12012 * costheta**6 + 6930 * costheta**4 -
            1260 * costheta**2 + 35)
    elif order == 1:
        return -0.1090412458987799555260481189135496447649 * sintheta * (
            715 * costheta**7 - 1001 * costheta**5 + 385 * costheta**3 -
            35 * costheta) * np.cos(phi)
    elif order == 2:
        return 0.4561522584349094705715788531131085781609 * sintheta**2 * (
            143 * costheta**6 - 143 * costheta**4 + 33 * costheta**2 -
            1) * np.cos(2 * phi)
    elif order == 3:
        return -1.2352661552955440759618898084529796356280 * sintheta**3 * (
            39 * costheta**5 - 26 * costheta**3 + 3 * costheta) * np.cos(
                3 * phi)
    elif order == 4:
        return 0.4784165247593306972508126646608462146527 * sintheta**4 * (
            65 * costheta**4 - 26 * costheta**2 + 1) * np.cos(4 * phi)
    elif order == 5:
        return -3.4499106220981080175178857353706008701770 * sintheta**5 * (
            5 * costheta**3 - costheta) * np.cos(5 * phi)
    elif order == 6:
        return 0.5323327660595425743111699467876298815112 * sintheta**6 * (
            15 * costheta**2 - 1) * np.cos(6 * phi)
    elif order == 7:
        return -2.9157066406993194754817395422284133340830 * sintheta**7 * costheta * np.cos(
            7 * phi)
    elif order == 8:
        return 0.7289266601748298688704348855571033335207 * sintheta**8 * np.cos(
            8 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def _spherical_harmonics_band10(order: int, costheta: np.ndarray,
                                sintheta: np.ndarray,
                                phi: np.ndarray) -> np.float32:
    if order == -10:
        return 0.7673951182219900125615933844975307826007 * sintheta**10 * np.sin(
            10 * phi)
    elif order == -9:
        return -3.4318952998917144349297140639050697246138 * sintheta**9 * costheta * np.sin(
            9 * phi)
    elif order == -8:
        return 0.5567269327204184018890957936038479115945 * sintheta**8 * (
            19 * costheta**2 - 1) * np.sin(8 * phi)
    elif order == -7:
        return -1.3636969112298053697585581756292731165454 * sintheta**7 * (
            19 * costheta**3 - 3 * costheta) * np.sin(7 * phi)
    elif order == -6:
        return 0.1653725413626187638731375305623700315483 * sintheta**6 * (
            323 * costheta**4 - 102 * costheta**2 + 3) * np.sin(6 * phi)
    elif order == -5:
        return -0.2958273952789690039984672558593464519425 * sintheta**5 * (
            323 * costheta**5 - 170 * costheta**3 + 15 * costheta) * np.sin(
                5 * phi)
    elif order == -4:
        return 0.4677441816782421841904245560595109212672 * sintheta**4 * (
            323 * costheta**6 - 255 * costheta**4 + 45 * costheta**2 -
            1) * np.sin(4 * phi)
    elif order == -3:
        return -0.6614901654504750554925501222494801261932 * sintheta**3 * (
            323 * costheta**7 - 357 * costheta**5 + 105 * costheta**3 -
            7 * costheta) * np.sin(3 * phi)
    elif order == -2:
        return 0.0648644473400325399686187038310999771570 * sintheta**2 * (
            4199 * costheta**8 - 6188 * costheta**6 + 2730 * costheta**4 -
            364 * costheta**2 + 7) * np.sin(2 * phi)
    elif order == -1:
        return -0.0748990122652081836719047053056755319005 * sintheta * (
            4199 * costheta**9 - 7956 * costheta**7 + 4914 * costheta**5 -
            1092 * costheta**3 + 63 * costheta) * np.sin(phi)
    elif order == 0:
        return 0.0050496903767836039532405338020352903762 * (
            46189 * costheta**10 - 109395 * costheta**8 + 90090 * costheta**6 -
            30030 * costheta**4 + 3465 * costheta**2 - 63)
    elif order == 1:
        return -0.0748990122652081836719047053056755319005 * sintheta * (
            4199 * costheta**9 - 7956 * costheta**7 + 4914 * costheta**5 -
            1092 * costheta**3 + 63 * costheta) * np.cos(phi)
    elif order == 2:
        return 0.0648644473400325399686187038310999771570 * sintheta**2 * (
            4199 * costheta**8 - 6188 * costheta**6 + 2730 * costheta**4 -
            364 * costheta**2 + 7) * np.cos(2 * phi)
    elif order == 3:
        return -0.6614901654504750554925501222494801261932 * sintheta**3 * (
            323 * costheta**7 - 357 * costheta**5 + 105 * costheta**3 -
            7 * costheta) * np.cos(3 * phi)
    elif order == 4:
        return 0.4677441816782421841904245560595109212672 * sintheta**4 * (
            323 * costheta**6 - 255 * costheta**4 + 45 * costheta**2 -
            1) * np.cos(4 * phi)
    elif order == 5:
        return -0.2958273952789690039984672558593464519425 * sintheta**5 * (
            323 * costheta**5 - 170 * costheta**3 + 15 * costheta) * np.cos(
                5 * phi)
    elif order == 6:
        return 0.1653725413626187638731375305623700315483 * sintheta**6 * (
            323 * costheta**4 - 102 * costheta**2 + 3) * np.cos(6 * phi)
    elif order == 7:
        return -1.3636969112298053697585581756292731165454 * sintheta**7 * (
            19 * costheta**3 - 3 * costheta) * np.cos(7 * phi)
    elif order == 8:
        return 0.5567269327204184018890957936038479115945 * sintheta**8 * (
            19 * costheta**2 - 1) * np.cos(8 * phi)
    elif order == 9:
        return -3.4318952998917144349297140639050697246138 * sintheta**9 * costheta * np.cos(
            9 * phi)
    elif order == 10:
        return 0.7673951182219900125615933844975307826007 * sintheta**10 * np.cos(
            10 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def _spherical_harmonics_band12(order: int, costheta: np.ndarray,
                                sintheta: np.ndarray,
                                phi: np.ndarray) -> np.float32:
    if order == -12:
        return 0.5662666637421911709547090258698966891658 * sintheta**12 * np.sin(
            12 * phi)
    elif order == -11:
        return -2.7741287690330965092773166877539934982424 * sintheta**11 * costheta * np.sin(
            11 * phi)
    elif order == -10:
        return 0.4090229723318165722310088483734983242130 * sintheta**10 * (
            23 * costheta**2 - 1) * np.sin(10 * phi)
    elif order == -9:
        return -1.1076394452006765536029569648662831371250 * sintheta**9 * (
            23 * costheta**3 - 3 * costheta) * np.sin(9 * phi)
    elif order == -8:
        return 0.3625601143107851015353138806669183470796 * sintheta**8 * (
            161 * costheta**4 - 42 * costheta**2 + 1) * np.sin(8 * phi)
    elif order == -7:
        return -0.7251202286215702030706277613338366941591 * sintheta**7 * (
            161 * costheta**5 - 70 * costheta**3 + 5 * costheta) * np.sin(
                7 * phi)
    elif order == -6:
        return 0.0679137317817836802249008221855060202705 * sintheta**6 * (
            3059 * costheta**6 - 1995 * costheta**4 + 285 * costheta**2 -
            5) * np.sin(6 * phi)
    elif order == -5:
        return -0.7623297485540852851293571095780914553217 * sintheta**5 * (
            437 * costheta**7 - 399 * costheta**5 + 95 * costheta**3 -
            5 * costheta) * np.sin(5 * phi)
    elif order == -4:
        return 0.0653692366445350742673043783584385579335 * sintheta**4 * (
            7429 * costheta**8 - 9044 * costheta**6 + 3230 * costheta**4 -
            340 * costheta**2 + 5) * np.sin(4 * phi)
    elif order == -3:
        return -0.0871589821927134323564058378112514105780 * sintheta**3 * (
            7429 * costheta**9 - 11628 * costheta**7 + 5814 * costheta**5 -
            1020 * costheta**3 + 45 * costheta) * np.sin(3 * phi)
    elif order == -2:
        return 0.1067475164362366128085636211368571894395 * sintheta**2 * (
            7429 * costheta**10 - 14535 * costheta**8 + 9690 * costheta**6 -
            2550 * costheta**4 + 225 * costheta**2 - 3) * np.sin(2 * phi)
    elif order == -1:
        return -0.0172039200193992376329375079035743398052 * sintheta * (
            52003 * costheta**11 - 124355 * costheta**9 + 106590 * costheta**7 -
            39270 * costheta**5 + 5775 * costheta**3 -
            231 * costheta) * np.sin(phi)
    elif order == 0:
        return 0.0013774159754583893724318345985370424459 * (
            676039 * costheta**12 - 1939938 * costheta**10 +
            2078505 * costheta**8 - 1021020 * costheta**6 +
            225225 * costheta**4 - 18018 * costheta**2 + 231)
    elif order == 1:
        return -0.0172039200193992376329375079035743398052 * sintheta * (
            52003 * costheta**11 - 124355 * costheta**9 + 106590 * costheta**7 -
            39270 * costheta**5 + 5775 * costheta**3 -
            231 * costheta) * np.cos(phi)
    elif order == 2:
        return 0.1067475164362366128085636211368571894395 * sintheta**2 * (
            7429 * costheta**10 - 14535 * costheta**8 + 9690 * costheta**6 -
            2550 * costheta**4 + 225 * costheta**2 - 3) * np.cos(2 * phi)
    elif order == 3:
        return -0.0871589821927134323564058378112514105780 * sintheta**3 * (
            7429 * costheta**9 - 11628 * costheta**7 + 5814 * costheta**5 -
            1020 * costheta**3 + 45 * costheta) * np.cos(3 * phi)
    elif order == 4:
        return 0.0653692366445350742673043783584385579335 * sintheta**4 * (
            7429 * costheta**8 - 9044 * costheta**6 + 3230 * costheta**4 -
            340 * costheta**2 + 5) * np.cos(4 * phi)
    elif order == 5:
        return -0.7623297485540852851293571095780914553217 * sintheta**5 * (
            437 * costheta**7 - 399 * costheta**5 + 95 * costheta**3 -
            5 * costheta) * np.cos(5 * phi)
    elif order == 6:
        return 0.0679137317817836802249008221855060202705 * sintheta**6 * (
            3059 * costheta**6 - 1995 * costheta**4 + 285 * costheta**2 -
            5) * np.cos(6 * phi)
    elif order == 7:
        return -0.7251202286215702030706277613338366941591 * sintheta**7 * (
            161 * costheta**5 - 70 * costheta**3 + 5 * costheta) * np.cos(
                7 * phi)
    elif order == 8:
        return 0.3625601143107851015353138806669183470796 * sintheta**8 * (
            161 * costheta**4 - 42 * costheta**2 + 1) * np.cos(8 * phi)
    elif order == 9:
        return -1.1076394452006765536029569648662831371250 * sintheta**9 * (
            23 * costheta**3 - 3 * costheta) * np.cos(9 * phi)
    elif order == 10:
        return 0.4090229723318165722310088483734983242130 * sintheta**10 * (
            23 * costheta**2 - 1) * np.cos(10 * phi)
    elif order == 11:
        return -2.7741287690330965092773166877539934982424 * sintheta**11 * costheta * np.cos(
            11 * phi)
    elif order == 12:
        return 0.5662666637421911709547090258698966891658 * sintheta**12 * np.cos(
            12 * phi)
    else:
        raise Exception("Failed to get correct order!")


@njit(cache=True)
def spherical_harmonics(band: int, order: int, costheta: np.ndarray,
                        sintheta: np.ndarray, phi: np.ndarray) -> np.float32:
    if band == 0:
        return _spherical_harmonics_band0(order, costheta, sintheta, phi)
    if band == 2:
        return _spherical_harmonics_band2(order, costheta, sintheta, phi)
    if band == 4:
        return _spherical_harmonics_band4(order, costheta, sintheta, phi)
    if band == 6:
        return _spherical_harmonics_band6(order, costheta, sintheta, phi)
    if band == 8:
        return _spherical_harmonics_band8(order, costheta, sintheta, phi)
    if band == 10:
        return _spherical_harmonics_band10(order, costheta, sintheta, phi)
    if band == 12:
        return _spherical_harmonics_band12(order, costheta, sintheta, phi)
