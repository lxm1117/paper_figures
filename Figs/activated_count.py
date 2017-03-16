import numpy as np
import matplotlib.pylab as plt	


def cam_count(mol,state):
	# cam only
	if state=='cam':
		return mol[2]+mol[4]+mol[5]+mol[13]+mol[6]+mol[8]+mol[14]+mol[16]
	# cam-K
	elif state=='camK':
		return mol[18]+mol[20]+mol[21]+mol[29]+mol[22]+mol[24]+mol[30]+mol[32]
	# cam-Kp
	elif state=='camKp':
		return mol[50]+mol[52]+mol[53]+mol[61]+mol[54]+mol[56]+mol[62]+mol[64]
	elif state=='camNg':
		return mol[66]+mol[68]+mol[69]+mol[77]+mol[70]+mol[72]+mol[78]+mol[80]

def camca_count(mol,state):
	if state=='cam':
		return mol[2]+ 2*mol[4]+ mol[5]+ 2*mol[13]+ 2*mol[6]+ 3*mol[8]+ 3*mol[14]+ 4*mol[16]
	elif state=='camK':	# 16+
		return mol[18]+2*mol[20]+ mol[21]+ 2*mol[29]+ 2*mol[22]+ 3*mol[24]+ 3*mol[30]+ 4*mol[32]
	# cam-Kp 48+
	elif state=='camKp':
		return mol[50]+ 2*mol[52]+ mol[53]+ 2*mol[61]+ 2*mol[54]+ 3*mol[56]+ 3*mol[62]+ 4*mol[64]
	elif state=='camNg':	# 64+
		return mol[66]+ 2*mol[68]+ mol[69]+ 2*mol[77]+ 2*mol[70]+ 3*mol[72]+ 3*mol[78]+ 4*mol[80]

def	camkii_count(mol,state):
# subunits that are only bound" 'b'-'bp'
# subunits that are autonomous 'p'-'bp'
	if state=='b':
		return mol[2]+mol[4]+mol[6]+mol[8]+mol[16]+mol[10]+mol[12]+mol[14]
	elif state=='bp':
		return mol[16]+mol[10]+mol[12]+mol[14]  
	elif state=='p':
		return mol[16]+mol[9]+mol[11]+mol[13]+mol[15]+mol[10]+mol[12]+mol[14]	
	elif state=='n': # neighbor subunits are either b or p
		return np.sum(mol[1:],0)-mol[2]-mol[9]-mol[10]
	elif state=='n2':	# both its neighbor and itself are activated
		return mol[4]+mol[6]+mol[8]+mol[16]+mol[11]+mol[13]+mol[15]+mol[12]+mol[14]

	elif state=='n2p': # both subunits are phosphorylated
		return mol[16]+mol[13]+mol[15]+mol[14]


def graph(ca1,cam1,bK1,ca2,cam2,bK2):
	fig=plt.figure(1)
	total_cam=np.sum(cam1[1:],0)
	total_bK=np.sum(bK1[1:],0)

	ax1=fig.add_subplot(431);
	ax1.plot(ca1[0],ca1[1],'g',ca2[0],ca2[1],'orange')
	ax1.set_ylabel('ca[1]')

	ax2=fig.add_subplot(432);
	ax2.plot(ca1[0],ca1[2],'g',ca2[0],ca2[2],'orange')
	ax2.set_ylabel('ca[2]')

	ax3=fig.add_subplot(433);
	ax3.plot(cam1[0],cam1[1]/total_cam,'g',cam2[0],cam2[1]/total_cam,'orange')
	ax3.set_ylabel('cam[1]')

	ax4=fig.add_subplot(434);
	ax4.plot(cam1[0],cam1[17]/total_cam,'g',cam2[0],cam2[17]/total_cam,'orange')
	ax4.set_ylabel('cam[17]')

	ax5=fig.add_subplot(435);
	ax5.plot(cam1[0],cam1[49]/total_cam,'g',cam2[0],cam2[49]/total_cam,'orange')
	ax5.set_ylabel('cam[49]')

	ax6=fig.add_subplot(436);
	ax6.plot(cam1[0],cam1[65]/total_cam,'g',cam2[0],cam2[65]/total_cam,'orange')
	ax6.set_ylabel('cam[65]')

	ax7=fig.add_subplot(437);
	ax7.plot(cam1[0],cam_count(cam1,'cam')/total_cam,'g',cam2[0],cam_count(cam2,'cam')/total_cam,'orange')
	ax7.set_ylabel('cam')

	ax8=fig.add_subplot(438);
	ax8.plot(cam1[0],cam_count(cam1,'camK')/total_cam,'g',cam2[0],cam_count(cam2,'camK')/total_cam,'orange')
	ax8.set_ylabel('camK')

	ax9=fig.add_subplot(439);
	ax9.plot(cam1[0],cam_count(cam1,'camKp')/total_cam,'g',cam2[0],cam_count(cam2,'camKp')/total_cam,'orange')
	ax9.set_ylabel('camKp')

	ax10=fig.add_subplot(4,3,10);
	ax10.plot(cam1[0],cam_count(cam1,'camNg')/total_cam,'g',cam2[0],cam_count(cam2,'camNg')/total_cam,'orange')
	ax10.set_ylabel('camNg')

	ax11=fig.add_subplot(4,3,11);
	ax11.plot(bK1[0],camkii_count(bK1,'n2'),'g',bK2[0],camkii_count(bK2,'n2'),'orange')
	ax11.set_ylabel('bK n2')

	ax12=fig.add_subplot(4,3,12);
	ax12.plot(bK1[0],camkii_count(bK1,'p'),'g',bK2[0],camkii_count(bK2,'p'),'orange')
	ax12.set_ylabel('bK p')
	
	plt.show()


def ca_graph(ca_n,ca_c):
	fig=plt.figure(1)
	ax1=fig.add_subplot(211)
	ax1.plot(ca_n[0],ca_n[2],ca_c[0],ca_c[2],ca_c[0],ca_c[2]+ca_n[2],'k--',lw=0.5)
	ax2=fig.add_subplot(212)
	ax2.plot(ca_n[0],ca_n[1],ca_c[0],ca_c[1],ca_c[0],ca_c[1]+ca_n[1],'k--',lw=0.5)
	#plt.show()

def cam_graph(cam1,cam2,flag):
	fig=plt.figure(1)
	ax=fig.add_subplot(111)
	if flag=='cam':
		ax.plot(cam1[0],cam_count(cam1,'cam'),'r',lw=0.5)
		ax.plot(cam1[0],cam_count(cam1,'camK'),'b',lw=0.5)
		ax.plot(cam1[0],cam_count(cam1,'camKp'),'g',lw=0.5)
	
		ax.plot(cam2[0],cam_count(cam2,'cam'),'r',lw=0.3)
		ax.plot(cam2[0],cam_count(cam2,'camK'),'b',lw=0.3)
		ax.plot(cam2[0],cam_count(cam2,'camKp'),'g',lw=0.3)
	elif flag=='ca':
		ax.plot(cam1[0],camca_count(cam1,'cam'),'r',lw=0.5)
		ax.plot(cam1[0],camca_count(cam1,'camK'),'b',lw=0.5)
		ax.plot(cam1[0],camca_count(cam1,'camKp'),'g',lw=0.5)
	
		ax.plot(cam2[0],camca_count(cam2,'cam'),'r',lw=0.3)
		ax.plot(cam2[0],camca_count(cam2,'camK'),'b',lw=0.3)
		ax.plot(cam2[0],camca_count(cam2,'camKp'),'g',lw=0.3)


def bK_graph(bK1,bK2,flag):
	fig=plt.figure(1)

	if flag=='b':
		f=-8
	elif flag=='p':
		f=0

	ax1=fig.add_subplot(241)
	ax1.plot(bK1[0],bK1[9+f],'r',bK2[0],bK2[9+f],'b')
	ax1.set_ylabel('0001')

	ax2=fig.add_subplot(242)
	ax2.plot(bK1[0],bK1[10+f],'r',bK2[0],bK2[10+f],'b')
	ax2.set_ylabel('1001')
	
	ax3=fig.add_subplot(243)
	ax3.plot(bK1[0],bK1[11+f],'r',bK2[0],bK2[11+f],'b')
	ax3.set_ylabel('0101')
	
	ax4=fig.add_subplot(244)
	ax4.plot(bK1[0],bK1[12+f],'r',bK2[0],bK2[12+f],'b')
	ax4.set_ylabel('1101')

	ax5=fig.add_subplot(245)
	ax5.plot(bK1[0],bK1[13+f],'r',bK2[0],bK2[13+f],'b')
	ax5.set_ylabel('0011')

	ax6=fig.add_subplot(246)
	ax6.plot(bK1[0],bK1[14+f],'r',bK2[0],bK2[14+f],'b')
	ax6.set_ylabel('1011')

	ax7=fig.add_subplot(247)
	ax7.plot(bK1[0],bK1[15+f],'r',bK2[0],bK2[15+f],'b')
	ax7.set_ylabel('0111')

	ax8=fig.add_subplot(248)
	ax8.plot(bK1[0],bK1[16+f],'r',bK2[0],bK2[16+f],'b')
	ax8.set_ylabel('1111')

	plt.show()
	return fig
