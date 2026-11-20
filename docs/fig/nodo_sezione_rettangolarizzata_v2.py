import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon
import numpy as np

BLUE='#2F6FD6'; BFILL='#DCE8FA'; OR='#E8692B'; OFILL='#FBE3D6'; GREY='#9A9A9A'; DK='#222222'; MID='#666666'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'mathtext.fontset':'dejavusans'})

fig=plt.figure(figsize=(7.2,10.2))
fig.text(0.03,0.985,"Il nodo si risolve in bankfull nella forma BRT (6 equazioni), poi si torna al volo",fontsize=10.5,weight='bold',va='top')
fig.text(0.03,0.962,"Versione aggiornata (20/11/2026): limiti del tirante bankfull, $S_0$ da finestra, pendenza di equilibrio $S_b$ incognita,\nritorno al volo con il moto uniforme",fontsize=7.6,color=MID,va='top')

def harrow(ax,x0,x1,y,txt,c=MID,dy=0.0,fs=8):
    ax.annotate('',xy=(x1,y),xytext=(x0,y),arrowprops=dict(arrowstyle='<|-|>',color=c,lw=0.8,mutation_scale=7))
    ax.text((x0+x1)/2,y+dy,txt,ha='center',va='bottom' if dy>=0 else 'top',fontsize=fs,color=DK)
def varrow(ax,x,y0,y1,c=MID):
    ax.annotate('',xy=(x,y1),xytext=(x,y0),arrowprops=dict(arrowstyle='<|-|>',color=c,lw=0.9,mutation_scale=7))

# ---------------- panel a ----------------
ax=fig.add_axes([0.03,0.655,0.94,0.27]); ax.set_xlim(0,100); ax.set_ylim(-14,62); ax.axis('off')
ax.text(0,61,"a) Madre: tirante bankfull limitato dalle sponde, rettangolo equivalente",weight='bold',fontsize=9.5,va='top')
px=[0,8,16,22,27,33,40,60,68,74,80,86,92,100]
pz=[44,46,44,37,30,14,8,8,14,26,38,48,50,50]
ax.plot(px,pz,color=GREY,lw=1.2)
ax.text(0,38,"profilo del DTM corretto\n(transetto ortogonale\nalla centerline)",color=GREY,fontsize=7.0,va='top')
eta0=8; DS=40; spill=44; Dc=18
# spill level = min of the side maxima (left max 46, right 50 -> 46)? set spill at 46 - 0.01
spill=46-0.6
D0=min(DS,spill-eta0); D0=max(D0,Dc)
# Shields level (above spill)
ax.plot([20,96],[eta0+DS,eta0+DS],ls=(0,(4,3)),color=BLUE,lw=0.9)
ax.text(62,eta0+DS+0.6,"$\\eta_0+D_S$ (Shields)",fontsize=7.3,color=BLUE,va='bottom')
ax.plot([4,96],[spill,spill],ls=':',color=DK,lw=0.9)
ax.text(87,spill-0.8,"$z_{sfioro}$",fontsize=7.5,va='top')
ax.plot([28,72],[eta0+Dc,eta0+Dc],ls=':',color=MID,lw=0.8)
ax.text(73,eta0+Dc,"$\\eta_0+D_c$",fontsize=7.5,color=MID,va='center')
# equivalent rectangle at level eta0+D0
lvl=eta0+D0
# bank intersections at level lvl: left between x=8..16 (46..44)... approximate
xL=17.2; xR=84.5
W=xR-xL
ax.add_patch(Rectangle((xL,eta0),W,D0,facecolor=BFILL,edgecolor=BLUE,lw=1.4,alpha=0.85,zorder=1))
ax.plot([xL,xR],[eta0,eta0],color=BLUE,lw=1.4)
# flight water level
wl=24
ax.plot([30.6,71.2],[wl,wl],color=OR,lw=1.6,zorder=3)
ax.text(34,wl-1.2,"pelo libero del volo ($z_{DTM}$)",color=OR,fontsize=7.3,va='top')
varrow(ax,52,eta0,lvl,BLUE)
ax.text(53,lvl-6,"$D_0$",fontsize=9,weight='bold')
ax.text(0,eta0+2,"$\\eta_0$ = fondo\nidraulico",fontsize=7.3,va='center')
harrow(ax,xL,xR,52.5,"$W_0^{bf}=A_0/D_0\\ \\geq\\ W_0^{m}$",BLUE)
ax.plot([xL,xL],[lvl,52.5],color=GREY,lw=0.5,ls=':'); ax.plot([xR,xR],[lvl,52.5],color=GREY,lw=0.5,ls=':')
harrow(ax,30.6,71.2,4.2,"$W_0^m$: larghezza bagnata al volo (maschera)",MID,dy=-1.2,fs=7.5)
ax.text(0,-3.5,"$D_0=\\max\\left[D_c,\\ \\min\\left(D_S,\\ z_{sfioro}-\\eta_0\\right)\\right]$,   $D_S=1.62\\,\\theta^*_c\\,\\Delta d_{50}/S_0$,   $D_c=\\theta^*_c\\,\\Delta d_{50}/S_0$",fontsize=8.2,va='top')
ax.text(0,-8.5,"$S_0$: retta sul pelo libero entro $L$ = 400 m a monte del nodo. $z_{sfioro}$: il più basso dei due punti più alti\ndel terreno tra $W^m/2$ e $3W^m$ dalla centerline, uno per lato.",fontsize=7.4,va='top',color=MID)

# ---------------- panel b ----------------
ax=fig.add_axes([0.03,0.355,0.94,0.29]); ax.set_xlim(0,100); ax.set_ylim(-22,62); ax.axis('off')
ax.text(0,61,"b) Nodo in bankfull (BRT): pelo libero comune, stessa pendenza $S_b$ nei rami (incognita)",weight='bold',fontsize=9.5,va='top')
# planform inset on the right
ix0=66
ax.add_patch(Rectangle((ix0,30),12,16,facecolor='#F3F6FB',edgecolor=BLUE,lw=0.9))
ax.plot([ix0,ix0+12],[38,38],ls=(0,(3,2)),color=BLUE,lw=0.7)
ax.text(ix0+8,42,"cella 1",ha='center',fontsize=6.8,color=DK); ax.text(ix0+8,33.2,"cella 2",ha='center',fontsize=6.8,color=DK)
ax.add_patch(Polygon([[ix0+12,46],[ix0+26,52],[ix0+26,44],[ix0+12,38]],closed=True,facecolor='#F3F6FB',edgecolor=BLUE,lw=0.9))
ax.add_patch(Polygon([[ix0+12,38],[ix0+26,30],[ix0+26,24],[ix0+12,30]],closed=True,facecolor='#F3F6FB',edgecolor=BLUE,lw=0.9))
ax.text(ix0+27,48,"ramo 1",fontsize=6.8,va='center'); ax.text(ix0+27,27,"ramo 2",fontsize=6.8,va='center')
ax.annotate('',xy=(ix0+6.5,35.5),xytext=(ix0+6.5,40.5),arrowprops=dict(arrowstyle='-|>',color=OR,lw=1,mutation_scale=7))
ax.text(ix0+1,37.6,"$q_{sy,x}$",fontsize=7,color=OR,va='center')
harrow(ax,ix0,ix0+12,26.5,"$\\alpha W_0$",MID,dy=-1,fs=7)
ax.text(ix0-1,53.5,"pianta: due celle a monte del nodo",fontsize=6.8,color=MID)
# section view
ws=70/100*60  # scale
x0=2; W1=33; W2=22
wsurf=50
D1=34; D2=22
eta1=wsurf-D1; eta2=wsurf-D2; eta0=wsurf-27
ax.add_patch(Rectangle((x0,eta1),W1,D1,facecolor=BFILL,edgecolor=BLUE,lw=1.4))
ax.add_patch(Rectangle((x0+W1,eta2),W2,D2,facecolor=BFILL,edgecolor=BLUE,lw=1.4))
ax.plot([x0,x0+W1+W2+4],[wsurf,wsurf],color=BLUE,lw=1.2)
ax.text(x0+W1+W2+4.5,wsurf,"pelo libero comune",fontsize=7.5,va='center',color=BLUE)
ax.text(x0+(W1+W2)/2,wsurf+3.5,"$\\eta_0+D_0=\\eta_1+D_1=\\eta_2+D_2$   →   $\\eta_1-\\eta_2=D_2-D_1$",ha='center',fontsize=8)
ax.plot([x0,x0+W1+W2+4],[eta0,eta0],ls=(0,(4,3)),color=GREY,lw=0.8)
ax.text(x0+W1+W2+4.5,eta0,"$\\eta_0$",fontsize=7.5,color=MID,va='center')
varrow(ax,x0+W1*0.75,eta1,wsurf-6,BLUE); ax.text(x0+W1*0.75+1,eta1+8,"$D_1$",fontsize=9,weight='bold'); ax.text(x0+W1*0.75+1,eta1+4.5,"incognita",fontsize=7,color=MID)
varrow(ax,x0+W1+W2*0.75,eta2,wsurf-6,BLUE); ax.text(x0+W1+W2*0.75-8,eta2+8,"$D_2$",fontsize=9,weight='bold'); ax.text(x0+W1+W2*0.75-8,eta2+4.5,"incognita",fontsize=7,color=MID)
ax.text(x0+1,wsurf-4,"$\\psi Q_0^{bf}=W_1 D_1^{5/3}\\sqrt{S_b}/n$",fontsize=7.6)
ax.text(x0+W1+0.8,wsurf-4,"$(1-\\psi)Q_0^{bf}=W_2 D_2^{5/3}\\sqrt{S_b}/n$",fontsize=7.0)
ax.text(x0+W1-1,eta1+1,"$\\eta_1$",fontsize=7.5,ha='right'); ax.text(x0+W1+W2-1,eta2+1,"$\\eta_2$",fontsize=7.5,ha='right')
harrow(ax,x0,x0+W1,eta1-4,"$W_1=W_0\\,W_1^{dtm}/(W_1^{dtm}+W_2^{dtm})$",MID,dy=-1.2,fs=7.2)
harrow(ax,x0+W1,x0+W1+W2,eta1-4,"$W_2=W_0-W_1$",MID,dy=-1.2,fs=7.2)
ax.text(0,4,"Incognite: $\\psi=Q_1/Q_0$, $D_1$, $D_2$, $S_b$, $f_{1F}$, $f_{2F}$. Equazioni: moto uniforme nei due rami (2),\nbilanci di cella per le classi F e C (2), continuità solida per F e C (2).",fontsize=7.6,va='top')
ax.text(0,-6.5,"$W_i^{dtm}$: larghezza bankfull alla prima sezione del ramo, misurata al livello di Shields calcolato con $S_0$.",fontsize=7.6,va='top')
ax.text(0,-11.5,"Non si usano più la chiusura $\\eta_0=(\\eta_1+\\eta_2)/2$ né le pendenze misurate dei rami. Scelta della soluzione:\nla bilanciata se è l'unica, altrimenti l'asimmetrica più vicina a $\\psi_{hyd}$.",fontsize=7.4,va='top',color=MID)

# ---------------- panel c ----------------
ax=fig.add_axes([0.03,0.01,0.94,0.3]); ax.set_xlim(0,100); ax.set_ylim(-8,62); ax.axis('off')
ax.text(0,61,"c) Ritorno al volo: portate ripartite con $\\psi$ e $\\psi_x$, fondo dei rami dal moto uniforme",weight='bold',fontsize=9.5,va='top')
x0=2; W1m=30; W2m=20
zt1=44; zt2=41
h1=17; h2=12
e1=zt1-h1; e2=zt2-h2
ax.plot([x0,x0+W1m+W2m+4],[52,52],ls=(0,(4,3)),color=GREY,lw=0.8)
ax.text(x0+W1m+W2m+4.5,52,"livello bankfull del nodo (riferimento)",fontsize=7.2,color=MID,va='center')
ax.add_patch(Rectangle((x0,e1),W1m,h1,facecolor=OFILL,edgecolor=OR,lw=1.4))
ax.add_patch(Rectangle((x0+W1m+2,e2),W2m,h2,facecolor=OFILL,edgecolor=OR,lw=1.4))
ax.plot([x0,x0+W1m],[zt1,zt1],color=OR,lw=1.8); ax.plot([x0+W1m+2,x0+W1m+2+W2m],[zt2,zt2],color=OR,lw=1.8)
ax.text(x0+W1m+W2m+5,zt1+1,"pelo libero del volo $z_{DTM,i}$\n(prima sezione di ciascun ramo)",fontsize=7.2,color=OR,va='center')
varrow(ax,x0+W1m-4,e1,zt1,OR); ax.text(x0+W1m-3.3,e1+h1/2,"$h_1$",fontsize=9,weight='bold',va='center')
varrow(ax,x0+W1m+2+W2m-4,e2,zt2,OR); ax.text(x0+W1m+2+W2m-3.3,e2+h2/2,"$h_2$",fontsize=9,weight='bold',va='center')
ax.text(x0+1,zt1-3.2,"$Q_1=\\psi\\,Q_0$",fontsize=8); ax.text(x0+1,zt1-7,"$Q_{s1,x}=\\psi_x\\,Q_{s0,x}$",fontsize=8)
ax.text(x0+W1m+3,zt2-3.2,"$Q_2=(1-\\psi)\\,Q_0$",fontsize=7.6); ax.text(x0+W1m+3,zt2-7,"$Q_{s2,x}=Q_{s0,x}-Q_{s1,x}$",fontsize=7.0)
ax.text(x0+1,e1+1,"$\\eta_1$",fontsize=7.5); ax.text(x0+W1m+3,e2+1,"$\\eta_2$",fontsize=7.5)
harrow(ax,x0,x0+W1m,e1-3.5,"$W_1^m$ (maschera)",MID,dy=-1.2,fs=7.5)
harrow(ax,x0+W1m+2,x0+W1m+2+W2m,e2-3.5,"$W_2^m$ (maschera)",MID,dy=-1.2,fs=7.5)
ax.text(x0+W1m+W2m+6,27,"$h_i=\\left(\\dfrac{Q_i\\,n}{W_i^m\\sqrt{S_i}}\\right)^{3/5}$",fontsize=9.5,va='center')
ax.text(x0+W1m+W2m+6,16,"$\\eta_i=z_{DTM,i}-h_i$",fontsize=9.5,va='center')
ax.text(0,12,"Ogni sezione dei rami, la prima compresa, si calcola come nel canale singolo:\nportata $Q_i$, larghezza della maschera $W_i^m$, pendenza del tratto $S_i$.",fontsize=7.6,va='top')
ax.text(0,1,"Le quote bankfull $\\eta_i=\\eta_0+D_0-D_i$ del pannello b e le larghezze riscalate $W_i\\,W_0^m/W_0^{bf}$ restano solo diagnostica:",fontsize=7.4,va='top',color=MID)
ax.text(0,-3,"usate al volo davano fondi sopra il pelo libero in 5–7 nodi.",fontsize=7.4,va='top',color=MID)
fig.savefig('/home/claude/fig/nodo_sezione_rettangolarizzata_v2.png',dpi=200,facecolor='white')
