from Core import UCDPhot


test = UCDPhot(path2data=r'G:\Kislovodsk\T1b_J-20200912\done')
test.aperture_photometry()
test.diff_photometry()


# files = glob.glob(r'G:\Kislovodsk\T1b_J-20200912\*.fts')
# with tqdm(total=len(files), desc='Astrometry.net on fits') as bar:
#     for item in files:
#         print('\n')
#         try:
#             test.do_astrometry(path2image=item, isFITS=0)
#         except Exception as e:
#             print(e.args)
#         bar.update(1)
#     bar.close()

# import numpy as np
# import matplotlib.pyplot as plt
#
#
# t = np.arange(0, 8.1, 0.1)
# plt.grid()
# vx = 10
# vy = 20
# plt.plot(vx*t, vy*t-5*t**2, label=f'v_x = {vx} м/с, v_y = {vy} м/с')
# # plt.plot(6*t, 40*t-5*t**2, label='v_x = 6 м/с, v_y = 40 м/с')
# # plt.plot(4*t, 45*t-5*t**2, label='v_x = 4 м/с, v_y = 45 м/с')
# plt.title('Траектория полёта камня в поле гравитации Земли, t = 8 с')
# plt.xlabel('x, м')
# plt.ylabel('y, м')
# plt.legend()
# plt.show()
