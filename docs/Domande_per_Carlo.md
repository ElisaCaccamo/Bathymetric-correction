# Domande per Carlo

## 19/11/2026

Per ognuna: il contesto, la domanda, e dove serve le opzioni considerate.

**1. Classi delle stazioni GNSS per la validazione.**

*Contesto.* Una delle sotto-domande del paper è dove il modulo morfodinamico migliora quello idraulico: nel canale singolo, a valle delle biforcazioni o a valle delle confluenze. Per rispondere, ognuna delle 48 stazioni GNSS (sezioni trasversali misurate) viene assegnata a un tipo di tratto, ricavato dalla topologia della rete. Con la rete corretta le stazioni sono: 37 nel canale singolo, 5 a valle di una biforcazione, 1 a valle di una confluenza, 5 in altri tratti a canali multipli. Nel canale singolo i due moduli hanno la stessa portata e lo stesso fondo, quindi lo stesso errore per costruzione: il confronto tra i moduli si gioca solo sulle altre 11 stazioni, e con 5, 1 e 5 stazioni per classe gli intervalli di confidenza sono molto larghi (per esempio RMSE 0.09–0.60 m negli altri tratti a canali multipli).

*Domanda.* Come conviene raggruppare le stazioni per il confronto tra i moduli?
- (a) Due classi: canale singolo (37) contro canali multipli (11).
- (b) Nessuna classe: differenza tra i moduli in funzione della distanza dal nodo a monte (11 stazioni).
- (c) Le quattro classi attuali, cioè canale singolo (37 stazioni), a valle di una biforcazione (5), a valle di una confluenza (1), altri tratti a canali multipli (5), riportate solo come descrizione: RMSE ed errore medio per classe, senza intervalli di confidenza né test statistici, perché con 1–5 stazioni non sarebbero significativi.

**2. Coefficiente $\theta^*_{BF}=1.62\,\theta^*_c$ e sponde basse.**

*Contesto.* Il tirante bankfull della madre è il tirante di Shields $D_S=1.62\,\theta^*_c\Delta d_{50}/S_0$, limitato dal livello di sfioro del transetto e mai sotto il tirante critico $D_c=\theta^*_c\Delta d_{50}/S_0$ (§9.1). Nel run del 18/11 sera il limite che decide è: tirante critico (sponde più basse di $D_c$) in 5 madri su 10, Shields pieno in 2, livello di sfioro in 3. La fonte del coefficiente 1.62 non è documentata; un riferimento possibile è Parker et al. (2007), che esprimono la geometria bankfull dei fiumi ghiaiosi a canale singolo tramite il rapporto tra Shields bankfull e Shields critico.

*Domanda.* Da quale riferimento viene 1.62? Vale per un fiume ghiaioso a canali multipli come l'Orco, o conviene calibrare $\theta^*_{BF}$ sulle sezioni a canale singolo con sponde ben definite? Limitare il tirante con il livello di sfioro è un criterio accettabile?

**3. Portata bankfull non coerente lungo la rete.**

*Contesto.* Ogni nodo calcola la sua portata bankfull localmente, con il moto uniforme della madre a larghezza, tirante e pendenza bankfull. I valori non sono coerenti lungo la rete: 15: 186, 106: 149, 205: 306, 244: 549, 469: 72, 584: 214, 645: 476, 1070: 23, 1293: 264, 1369: 1118 m³/s (a 1369 il tirante bankfull è 4.0 m).

*Domanda.* Conviene imporre una sola portata formativa all'ingresso della rete e ripartirla con $\psi$ a ogni nodo, così che la portata bankfull sia coerente lungo la rete? Oppure si accettano valori locali, dichiarando che ogni nodo è trattato come isolato? Se si impone una portata unica, con quale criterio sceglierla alla stazione di San Benigno?

**4. Regola di scelta tra le soluzioni asimmetriche.**

*Contesto.* Quando il nodo ha soluzioni asimmetriche, si adotta quella con $\psi$ più vicina a $\psi_{hyd}$, la quota di portata del ramo 1 nel modulo idraulico al volo (§9.1). Nei nodi 205, 584 e 645 nessuna soluzione è vicina: $\psi$ = 0.47 contro $\psi_{hyd}$ = 0.74; 0.51 contro 0.66; 0.46 contro 0.35.

*Domanda.* $\psi_{hyd}$, che descrive la ripartizione a portata ordinaria, è un riferimento adeguato per scegliere tra equilibri che si formano a bankfull? Alternative possibili: scegliere la soluzione in cui il ramo più basso nel DTM è quello che porta più portata, oppure quella coerente con il ramo più largo nella maschera.

**5. Ruolo del trasporto solido alla portata del volo.**

*Contesto.* Il trasporto solido entra nel modulo morfodinamico alle condizioni bankfull: al nodo, il sistema BRT (moto uniforme nei rami, bilanci di cella e continuità solida per le classi fine e grossolana) fissa la ripartizione della portata $\psi$, i tiranti e la composizione dei rami. Poi si torna alla portata del volo (10.8 m³/s) per ricostruire il fondo: la portata del volo si ripartisce con la stessa $\psi$ e il tirante di ogni sezione si calcola con il moto uniforme. A questa portata i flussi solidi calcolati sono dell'ordine di 10⁻¹⁶–10⁻²¹ m³/s, cioè nulli, molto sotto la soglia di moto. Il trasporto solido determina quindi il fondo al volo solo attraverso $\psi$; le ripartizioni dei flussi solidi al volo, $\psi_x$, agiscono su quantità nulle.

*Domanda.* È coerente con l'impostazione del metodo (configurazione del nodo formata a bankfull, ritorno al volo solo per la ripartizione della portata)? Va presentato così nel paper?

**6. Densità dei punti ADCP nella validazione.**

*Contesto.* Ai 48 transetti GNSS (rilievo a guado) si aggiungono 8 transetti ADCP dello stesso giorno (07/03/2025), che entrano nella validazione come stazioni come le altre, ognuna confrontata con la sua quota media. L'ADCP registra un punto ogni 0.1–0.3 m e resta fermo qualche secondo all'inizio dei transetti (49–271 punti per transetto); il GNSS ha un punto ogni 1.5 m circa (mediana 20 punti per stazione). Nelle statistiche calcolate su tutti i punti insieme (RMSE, errore medio) un transetto ADCP da 271 punti pesa quanto più di dieci stazioni GNSS.

*Aggiornamento dal run del 19/11/2026.* Gli ADCP sono 8 stazioni su 59 ma il 51% dei punti, e pesano 5–6 stazioni GNSS anche nelle statistiche di gruppo con il bootstrap delle stazioni. Il loro residuo è 0.3–0.6 m più alto di quello delle stazioni GNSS vicine: i transetti ADCP coprono soprattutto le parti profonde, dove non si poteva guadare, quindi la loro media di stazione è più bassa della media della sezione intera con cui si confronta il fondo corretto. Con gli ADCP l'errore medio su tutte le stazioni sale da +0.22 a +0.33 m.

*Domanda.* Conviene diradare i punti ADCP a una spaziatura simile a quella dei GNSS (per esempio 1.5 m), oppure dare lo stesso peso a ogni stazione nelle statistiche, oppure tenere tutti i punti? Visto che coprono solo la parte profonda, conviene confrontarli con la media di stazione come le sezioni GNSS, tenerli separati, oppure usarli solo dove completano una sezione GNSS?
