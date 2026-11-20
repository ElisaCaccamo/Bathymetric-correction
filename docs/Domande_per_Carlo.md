# Domande per Carlo

## 19/11/2026

Per ognuna: il contesto, la domanda, e dove serve le opzioni considerate.

**1. Classi delle stazioni GNSS per la validazione.**

*Contesto.* Una delle sotto-domande del paper è dove il modulo morfodinamico migliora quello idraulico: nel canale singolo, a valle delle biforcazioni o a valle delle confluenze. Per rispondere, ognuna delle 48 stazioni GNSS (sezioni trasversali misurate) viene assegnata a un tipo di tratto, ricavato dalla topologia della rete. Con la rete corretta le stazioni sono: 37 nel canale singolo, 5 a valle di una biforcazione, 1 a valle di una confluenza, 5 in altri tratti a canali multipli. Nel canale singolo i due moduli hanno la stessa portata e lo stesso fondo, quindi lo stesso errore per costruzione: il confronto tra i moduli si gioca solo sulle altre 11 stazioni, e con 5, 1 e 5 stazioni per classe gli intervalli di confidenza sono molto larghi (per esempio RMSE 0.09–0.60 m negli altri tratti a canali multipli).

*Domanda.* Come conviene raggruppare le stazioni per il confronto tra i moduli?
- (a) Due classi: canale singolo (37) contro canali multipli (11).
- (b) Nessuna classe: differenza tra i moduli in funzione della distanza dal nodo a monte (11 stazioni).
- (c) Le quattro classi attuali, cioè canale singolo (37 stazioni), a valle di una biforcazione (5), a valle di una confluenza (1), altri tratti a canali multipli (5), riportate solo come descrizione: RMSE ed errore medio per classe, senza intervalli di confidenza né test statistici, perché con 1–5 stazioni non sarebbero significativi.

**2. Coefficiente $`\theta^*_{BF}=1.62\,\theta^*_c`$ e sponde basse.**

*Contesto.* Il tirante bankfull della madre è il tirante di Shields $`D_S=1.62\,\theta^*_c\Delta d_{50}/S_0`$, limitato dal livello di sfioro del transetto e mai sotto il tirante critico $`D_c=\theta^*_c\Delta d_{50}/S_0`$ (§9.1). Nel run del 18/11 sera il limite che decide è: tirante critico (sponde più basse di $D_c$) in 5 madri su 10, Shields pieno in 2, livello di sfioro in 3. La fonte del coefficiente 1.62 non è documentata; un riferimento possibile è Parker et al. (2007), che esprimono la geometria bankfull dei fiumi ghiaiosi a canale singolo tramite il rapporto tra Shields bankfull e Shields critico.

*Domanda.* Da quale riferimento viene 1.62? Vale per un fiume ghiaioso a canali multipli come l'Orco, o conviene calibrare $`\theta^*_{BF}`$ sulle sezioni a canale singolo con sponde ben definite? Limitare il tirante con il livello di sfioro è un criterio accettabile?

**3. Portata bankfull non coerente lungo la rete.**

*Contesto.* Ogni nodo calcola la sua portata bankfull localmente, con il moto uniforme della madre a larghezza, tirante e pendenza bankfull. I valori non sono coerenti lungo la rete: 15: 186, 106: 149, 205: 306, 244: 549, 469: 72, 584: 214, 645: 476, 1070: 23, 1293: 264, 1369: 1118 m³/s (a 1369 il tirante bankfull è 4.0 m).

*Domanda.* Conviene imporre una sola portata formativa all'ingresso della rete e ripartirla con $\psi$ a ogni nodo, così che la portata bankfull sia coerente lungo la rete? Oppure si accettano valori locali, dichiarando che ogni nodo è trattato come isolato? Se si impone una portata unica, con quale criterio sceglierla alla stazione di San Benigno?

**4. Regola di scelta tra le soluzioni asimmetriche.**

*Contesto.* Quando il nodo ha soluzioni asimmetriche, si adotta quella con $\psi$ più vicina a $\psi_{hyd}$, la quota di portata del ramo 1 nel modulo idraulico al volo (§9.1). Nei nodi 205, 584 e 645 nessuna soluzione è vicina: $\psi$ = 0.47 contro $\psi_{hyd}$ = 0.74; 0.51 contro 0.66; 0.46 contro 0.35.

*Domanda.* $\psi_{hyd}$, che descrive la ripartizione a portata ordinaria, è un riferimento adeguato per scegliere tra equilibri che si formano a bankfull? Alternative possibili: scegliere la soluzione in cui il ramo più basso nel DTM è quello che porta più portata, oppure quella coerente con il ramo più largo nella maschera.

**5. Ruolo del trasporto solido alla portata del volo.**

*Contesto.* Il trasporto solido entra nel modulo morfodinamico alle condizioni bankfull: al nodo, il sistema BRT (moto uniforme nei rami, bilanci di cella e continuità solida per le classi fine e grossolana) fissa la ripartizione della portata $\psi$, i tiranti e la composizione dei rami. Poi si torna alla portata del volo (10.8 m³/s) per ricostruire il fondo: la portata del volo si ripartisce con la stessa $\psi$ e il tirante di ogni sezione si calcola con il moto uniforme. A questa portata i flussi solidi calcolati sono dell'ordine di 10⁻¹⁶–10⁻²¹ m³/s, cioè nulli, molto sotto la soglia di moto. Il trasporto solido determina quindi il fondo al volo solo attraverso $\psi$; le ripartizioni dei flussi solidi al volo, $\psi_x$, agiscono su quantità nulle.

*Domanda.* È coerente con l'impostazione del metodo (configurazione del nodo formata a bankfull, ritorno al volo solo per la ripartizione della portata)? Va presentato così nel paper?

**6. Dati ADCP: copertura dei transetti e peso nella validazione.**

*Contesto: come si fa la validazione.* Il modello ricostruisce in ogni sezione un fondo piatto, cioè una sola quota per tutta la larghezza bagnata, ottenuta come pelo libero meno il tirante medio della sezione. Per confrontarlo con i dati, i punti di ogni sezione misurata (stazione) vengono mediati, e ogni punto è confrontato con la media della sua stazione. Il confronto funziona se la stazione copre la sezione da sponda a sponda: allora la sua media è il fondo medio della sezione, la stessa grandezza che il modello ricostruisce. Le 48 sezioni GNSS rilevate a guado fanno così.

*Contesto: il problema degli ADCP.* Ai dati GNSS si aggiungono 8 transetti ADCP dello stesso giorno (07/03/2025). L'ADCP è stato usato dove non si poteva guadare, quindi i transetti coprono soprattutto la parte profonda della sezione e non arrivano alle sponde, dove il fondo risale. La media di un transetto ADCP è quindi più bassa del fondo medio della sezione: confrontata con il fondo piatto del modello, dà un residuo più alto anche se modello e misure sono corretti.

- Le quote ADCP sono affidabili: nei punti vicini a punti GNSS concordano entro 4–9 cm.
- Nel run del 19/11 il residuo medio delle stazioni ADCP è 0.3–0.6 m più alto di quello delle stazioni GNSS vicine (per esempio tratto 35: a06 +0.74 m contro s022 +0.05 m).
- Con gli ADCP l'errore medio su tutte le stazioni sale da +0.22 a +0.33 m e l'RMSE da 0.30 a 0.40 m: è un effetto della copertura, non dell'accuratezza del modello.

*Contesto: il peso dei punti.* L'ADCP registra un punto ogni 0.1–0.3 m e resta fermo qualche secondo all'inizio dei transetti (49–271 punti per transetto); il GNSS ha un punto ogni 1.5 m circa (mediana 20 punti per stazione). Gli 8 transetti sono il 14% delle stazioni ma il 51% dei punti, e le statistiche (RMSE, errore medio) sono calcolate sui punti: un transetto ADCP pesa quanto 5–6 stazioni GNSS, anche nel bootstrap delle stazioni.

*Domanda.* Come conviene usare gli ADCP nella validazione?
- (a) Tenerli insieme alle sezioni GNSS con la media di stazione, dichiarando l'effetto della copertura.
- (b) Tenerli separati: validazione GNSS per il fondo medio della sezione, e un confronto a parte per gli ADCP (per esempio punto per punto, o solo sulla quota minima del fondo).
- (c) Usarli solo dove completano una sezione GNSS vicina (la parte profonda che non si poteva guadare), per ottenere sezioni complete da sponda a sponda.

E per il peso: diradare i punti ADCP a una spaziatura simile a quella dei GNSS (circa 1.5 m), dare lo stesso peso a ogni stazione nelle statistiche, oppure tenere tutti i punti?

## 20/11/2026

**7. Volo 2019 e DTM corretto del 2019 nel paper.**

*Contesto.* Il lavoro usa il LiDAR 2025 (volo del 17/03/2025). Esiste anche un DTM del 2019 con la batimetria già corretta da altri autori. La correzione dovrebbe essere basata sulla stessa teoria del modulo idraulico (moto uniforme), ma i dettagli non sono ancora documentati: portata al volo, scabrezza, sezioni, ripartizione ai nodi, validazione. L'idea è usarlo solo per un breve paragrafo qualitativo: evoluzione dell'alveo tra il 2019 e il 2025 e dimostrazione delle potenzialità del metodo applicato a più anni. Non sarebbe una baseline.

*Domanda.* Il volo 2019 e il suo DTM corretto entrano nel paper? Se sì, come:
- (a) Solo il paragrafo qualitativo, con il DoD 2025 − 2019 e i limiti dichiarati (correzioni fatte da autori diversi).
- (b) Nessun confronto in questo paper: il multi-temporale resta uno sviluppo futuro, citato nelle conclusioni.

Se (a): chi può fornire i dettagli della correzione 2019?

**8. Piano dell'analisi di sensitività: parametri, intervalli e metriche.**

*Contesto.* Il run di riferimento è quello del 20/11 sera: 48 stazioni GNSS, RMSE 0.30 m, errore medio +0.22 m in entrambi i moduli. La proposta è di variare un parametro alla volta attorno ai valori di riferimento, tenendo fissi tutti gli altri, per circa 25 run in tutto. Per ogni run si guardano: RMSE ed errore medio di entrambi i moduli sull'insieme di stazioni della validazione principale (vedi domande 1 e 6); la differenza di RMSE tra modulo morfodinamico e idraulico; la quota di portata $\psi$ ai nodi e il numero di nodi asimmetrici; il P95 della mappa MM − HM.

| Parametro | Riferimento | Valori proposti | Motivo |
| --- | --- | --- | --- |
| $n$ di Manning | 0.030 | 0.025, 0.035, 0.040 | intervallo tipico dei fiumi ghiaiosi; il valore non è calibrato |
| Portata del volo | 10.76 m³/s | ±10% (9.7, 11.8) | incertezza della scala di deflusso |
| Passo tra le sezioni `dx` | 10 m | 5, 20 m | discretizzazione |
| Soglia dei salti isolati `glitch_tol` | 0.20 m | 0.10, 0.30 m | pre-processing di $z_{DTM}$ |
| $`\theta^*_{BF}/\theta^*_c`$ (`bankfull_factor`) | 1.62 | 1.2, 2.0 | fonte non documentata (domanda 2) |
| $`\theta^*_c`$ | 0.030 | 0.045 | soglia di moto |
| Coefficiente di Ikeda $r$ | 0.65 | 0.4, 0.9 | flusso laterale al nodo |
| Lunghezza delle celle $\alpha$ | 4 | 2, 6 | modello a due celle |
| Larghezza al nodo | bankfull dal DTM | maschera al volo | già provato il 18/11 (§9.3) |
| Finestra della pendenza al nodo `L_node` | 400 m | 200, 800 m | pendenza ai nodi |
| Soglia minima della pendenza al nodo `S_floor_prctile` | 0.05 | 0.10 | comprende il nodo di testa 13 (deciso il 20/11: soglia minima, dichiarata come limite) |
| Regola di scelta della soluzione | più vicina a $\psi_{hyd}$ | alternative della domanda 4 | equilibri multipli |

*Domanda.* Approvi parametri, intervalli e metriche? Ci sono parametri da aggiungere o da togliere? Basta variare un parametro alla volta, oppure servono combinazioni (per esempio $n$ insieme alla portata, che agiscono entrambe sul tirante)?

**9. Modifiche al metodo al nodo: da rivedere.**

*Contesto.* Rispetto alla formulazione iniziale (§4.3 del documento di teoria) il trattamento delle biforcazioni nel modulo morfodinamico è cambiato in cinque punti: tirante bankfull della madre limitato dalle sponde e mai sotto il tirante critico; pendenza della madre $S_0$ da una retta sul pelo libero a monte del nodo; tirante per le larghezze dei rami calcolato con $S_0$; sistema del nodo nella forma BRT a 6 equazioni, con pendenza di equilibrio $S_b$ incognita e senza la chiusura $\eta_0=(\eta_1+\eta_2)/2$; ritorno al volo con il fondo dei rami dal moto uniforme (larghezza della maschera, pendenza del tratto) invece che dalla quota $\eta_i$ del nodo. Le domande 2, 4 e 5 riguardano singoli aspetti di queste modifiche.

La §9.1 del documento di teoria ([Modifiche al metodo](TeoriaMetodologica_moduli_idraulico_morfodinamico.md#91-modifiche-al-metodo)) contiene la figura aggiornata del nodo (`fig/nodi_new.png`, da confrontare con `fig/nodi.png` del §4.3) e una tabella che riassume, per ciascuna modifica, la formulazione iniziale, il metodo attuale nel codice e il motivo del cambiamento.

*Domanda.* Le modifiche sono coerenti con l'impostazione del metodo e con BRT / Ragno et al. (2023)? C'è qualcosa che cambieresti prima dell'analisi di sensitività?
