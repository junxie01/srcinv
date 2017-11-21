from obspy.imaging.beachball import beachball

#(M11, M22, M33, M12, M13, M23)

mt_prior = [-0.7787, -4.279, 5.0577, -0.176, 1.465, 4.344]
beachball(mt_prior, size=200, linewidth=2, facecolor='b')

mt_post = [-0.56, -3.07, 12.48, 0.026, -4.86, 0.959]
beachball(mt_post, size=200, linewidth=2, facecolor='b')
