from Core import UCDPhot


bad_star_mask = [0]
test = UCDPhot(path2data=r'G:\Kislovodsk\T1c_J-20200913\done')
# test.do_astrometry(astrometry_net_API='hipfhzhlzygnlvix', bad_star_mask=bad_star_mask,
#                    path2images=r'G:\Kislovodsk\T1c_J-20201017')
test.aperture_photometry(cat_path=r'Cat.txt')
# test.create_humidity('2020-10-17.json')
# test.magic()
test.draw_curve()
# test.draw_magic_curve()
# test.plot_field()
