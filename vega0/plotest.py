import matplotlib.pyplot as plt 
from matplotlib import gridspec

plt.figure(figsize = (16,12))
gs1 = gridspec.GridSpec(2, 2, height_ratios=[2,.5])
gs1.update(wspace=0.0, hspace=0.0)

ax0 = plt.subplot(gs1[0])
ax0.plot_date(mwadates1, night1_ampscales, 'o', linestyle = ':', color='g')
ax0.errorbar(mwadates1, night1_ampscales, yerr=night1_stds, linestyle='None')
ax0.set_ylim([0.6, 1.4])
ax0.set_ylabel('Median amplitude scaling factor')
ax0.set_xlabel('Obsid date, time in 10/2013')
ax0.set_title('%s, %s Jy, mean qam: %s'%(source,str(round(np.nanmean(night1_flux),2)), str(round(night1_mdnmetric,2))))
ax0.grid()

ax1 = plt.subplot(gs1[1])#, sharey=ax0)
ax1.plot_date(mwadates2, night2_ampscales, 'o', linestyle = ':', color='g')
ax1.errorbar(mwadates2, night2_ampscales, yerr=night2_stds, linestyle='None')
ax1.set_ylim([0.6, 1.4])
ax1.set_ylabel('Median amplitude scaling factor')
ax1.set_xlabel('Obsid date, time in 10/2013')
ax2.set_title('%s, %s Jy, mean qam: %s'%(source,str(round(np.nanmean(night2_flux),2)), str(round(night2_mdnmetric,2))))
ax1.grid()






ax2 = plt.subplot(gs1[2])
ax2.plot(mwadates1, night1_err, color='g', linestyle='--')


ax3 = plt.subplot(gs1[3])
ax2.plot(mwadates2, night2_err, color='g', linestyle='--')
#plt.axis('on')

plt.setp(ax1.get_yticklabels(), visible=False)
plt.setp(ax3.get_yticklabels(), visible=False)

yticks = ax2.yaxis.get_major_ticks()
yticks[-1].label1.set_visible(False)

xticks = ax0.xaxis.get_major_ticks()
xticks[0].label1.set_visible(False)

xticks1 = ax1.xaxis.get_major_ticks()
xticks1[-1].label1.set_visible(False)



xticks2 = ax2.xaxis.get_major_ticks()
xticks2[-1].label1.set_visible(False)

ax0.spines['right'].set_visible(False)
ax1.spines['left'].set_visible(False)

ax2.spines['right'].set_visible(False)
ax3.spines['left'].set_visible(False)

ax1.yaxis.tick_left()
ax3.yaxis.tick_left()
                                                                                                                                                                        
plt.show()