# ClusteringSpartaco
di Sara Agavnì Castiglioni.


In questo pacchetto sono implementate la versione per il clustering, base e penalizzata, del modello di co-clustering Spartaco di A. Sottosanti e D. Risso.
Per maggiori dettagli si rimanda all'articolo del metodo spartaco https://arxiv.org/abs/2110.04872 e alla tesi di Laurea Magistrale di Sara Agavnì Castiglioni in Scienze Statistiche presso l'Università degli Studi di Padova. 


# Installazione
```
remotes::install_github("sacastiglioni/ClusteringSpartaco")
```

# Stima del modello
Sia ```x``` la matrice contenente l'espressione di ```nrow(x)``` geni misurata su ```ncol(x)``` spot. Le coordinate spaziali degli spot sono contenute nelle matrice ```coordinates```, mentre le etichette degli spot sono specificate nel vettore ```column.labels```. La funzione utilizza i parametri ```lambda.mu``` e ```lambda.tau``` per regolare la penalizzazione della stima del modello, che sono di default impostati a 0. È possibile eseguire la ricerca di ```K``` cluster di geni eseguendo il seguente codice:

```
library(ClusteringSpartaco)
RCspartaco(x = x, coordinates = coordinates, K = K, column.labels = column.labes)
```

# Convergenza
Per impostazione predefinita, la procedura di stima viene eseguita per un massimo di ```max.iter``` iterazioni, ma viene precedentemente interrotta se viene raggiunto un certo criterio di convergenza (```conv.criterion```). Se ```conv.criterion = NULL```, non c'è una condizione di fine e la procedura viene eseguita per ```max.iter``` iterazioni, altrimenti viene interrotta quando l'incremento della log-verosimiglianza di classificazione è inferiore a una certa soglia ```conv.criterion$epsilon``` per ```conv.criterion$iterations``` volte di seguito.
