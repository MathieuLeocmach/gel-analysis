import sys, os.path
from matplotlib.pyplot import *
from matplotlib.colors import LogNorm

#'2_ageing/broken_length_q2b.npy'
hq2blen = np.load(sys.argv[1])
fig = figure('broken length q2b')
clf()
imshow(
    hq2blen[3:,1:-1], 
    aspect='auto', interpolation='none',
    norm=LogNorm(vmin=1, vmax=hq2blen.max()), 
    extent=(np.linspace(0, 1,100)[1], np.linspace(0, 1,100)[-2], len(hq2blen)-1, 2)
    )
xlabel(r'$q_2$')
ylabel('Post-break distance')
colorbar()
#gca().annotate(
 #   "", 
  #  xy=(0.4, 12), xycoords='data', 
   # xytext=(np.linspace(0, 1,100)[1+np.argmax(hq2blen[3,1:-1])], 2.5), textcoords='data', 
    #arrowprops = dict(arrowstyle="fancy", color='1', connectionstyle="arc3,rad=-0.3"),
    #)

#longest length with at least 100 broken bonds
mml = np.where(hq2blen[3:,1:-1].sum(1)<100)[0][0]+3
step(np.linspace(0,1,100)[1:][[
    np.where(np.cumsum(l)>0.5*l.sum())[0][0] 
    for l in hq2blen[3:mml,1:-1]
    ]], np.arange(3,mml)-1, color='gray')
ylim(np.where(hq2blen[3:,1:-1].sum(1)<2)[0][0]+4,2)
draw()
savefig(os.path.splitext(sys.argv[1])[0]+'brokenlength_q2b.pdf')
