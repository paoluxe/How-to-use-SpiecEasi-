# Comment utiliser Spiec-Easi ?

---
titre: "Comment utiliser Spiec-Easi ?"   
auteure: "Paola Fournier"   
date: "21/06/2021"   
sortie: html_document   
---
# Comment utiliser SPIEC-EASI (Sparse and Compositionally Robust Inference of Microbial Ecological Networks)?  

Ce tutoriel s'adresse aux personnes novices voulant s'initier à l'inférence de réseau via le modèle Spiec-Easi (Kurtz et al., 2015).       

  Github : https://github.com/zdk123/SpiecEasi  
doi : https://doi.org/10.1371/journal.pcbi.1004226 

Mon choix s'est porté sur ce modèle car il s'attaque à deux problèmes majeurs des modèles basés sur la corrélation :   

  1) il utilise le concept de dépendance conditionnelle pour éviter la détection de taxons corrélés mais indirectement connectés.  La corrélation est une métrique par paire et donc limitée dans un cadre multivarié.  
  2) il tient compte de la compositionnalité des données en appliquant une transformation "clr" sur la table d'abondance d'entrée.

## **Quelques concepts importants avant de se lancer dans l'inférence de réseaux microbiens**  
### Compositionnalité des données   
Je voudrais d'abord parler de la façon dont les données de séquençage d'amplicons sont générées et pourquoi cela peut poser problème.     

**Qu'est-ce qu'une donnée compositionnelle?**   
es données compositionnelles sont définies comme un vecteur de nombres réels strictement positifs avec un total inconnu ou non informatif (par exemple la profondeur de séquençage) car l'abondance de chaque composant ne porte que des informations relatives (Pawlowsky-Glahn et al., 2015).   

**En quoi les données de séquençage d'amplicons sont compositionnelle?**      
Les instruments de séquençage ne peuvent fournir des lectures que jusqu'à leur capacité et, par conséquent, chaque échantillon est soumis à une contrainte arbitraire de somme constante. Les jeux de données dérivés des NGS sont donc de nature compositionnelle.  

**Quelles en sont les conséquences ?**   
Le nombre de lectures par échantillon n'est pas interprétable en soi : le nombre total de lectures attribuées à un échantillon ou à un taxa ne peut fournir aucune information sur le nombre de molécules dans l'échantillon d'origine.      
Nous pouvons seulement étudier les changements relatifs : nous ne pouvons pas comparer directement les échantillons parce qu'ils ont une profondeur de séquençage spécifique.     
De nombreuses méthodes statistiques ne doivent pas être utilisées, car elles supposent l'indépendance entre les caractéristiques, ce qui n'est pas le cas pour les données générées par le NGS.   

**Solutions proposées dans la littérature**   
Une technique populaire pour recalculer les abondances absolues est la comparaison du rapport logarithmique des comptages par rapport à une espèce de référence. Dans ce cas, l'espèce de référence est connue pour avoir une abondance approximativement stable à travers les populations. Ces rapports sont comparés au lieu des comptages directement.      
Si une telle espèce n'est pas connue ou n'est pas disponible, la référence peut être remplacée par une mesure composite robuste obtenue à partir de diverses espèces. Une de ces mesures est la moyenne géométrique de tous les comptages de l'échantillon. Cette méthode suppose qu'un agrégat de diverses espèces ne change pas en masse entre leurs environnements d'origine. 
![alt text](https://github.com/paoluxe/How-to-use-SpiecEasi-/blob/420285dcb6df1568a4d9084f0f5db615da09e6eb/Pictures%20README/alr.PNG)
![alt text](https://github.com/paoluxe/How-to-use-SpiecEasi-/blob/420285dcb6df1568a4d9084f0f5db615da09e6eb/Pictures%20README/clr.PNG)

### Principe des réseaux de co-occurrence    
Les réseaux de co-occurrence nécessitent comme données d’entrée des tables d’abondances de taxons provenant de multiples échantillons. Lorsque deux taxons co-occurrent ou montrent un schéma d’abondance similaire parmi plusieurs échantillons, une relation positive est supposée entre eux. Un lien positif est illustré en vert ici ou un 1 dans la table d'adjacence. Inversement lorsque ces deux taxons s’excluent mutuellement, une relation négative est assumée, illustré ici par un lien rouge reliant deux noeuds ou un -1 dans la table d'adjacence. 
![alt text](https://github.com/paoluxe/How-to-use-SpiecEasi-/blob/ef1049e3acef34d2c039fd7aee9abb8fbc162609/Pictures%20README/Principe%20R%C3%A9seaux%20de%20co-occurrence.PNG)

### Fonctionnement de SpiecEasi (méthode mb)  
Je décris ici la méthode "mb" i.e. la sélection de voisinage de Meinshausen-Buhlmann (Meinshausen et Buhlmann, 2006) qui constitue une des deux méthodes d'inférence de réseaux supportées par la fonction SpiecEasi. La méthode "mb" a montré de bons résultats et est rapide d'utilisation en comparaisonb avec la méthode "glasso" (Kurtz et al., 2015 ; Röttjers et Faust, 2018).   

A partir d’une table d’abondance de taxons sur de multiples échantillons, l'outil applique d'abord une transformation "clr" afin de recalculer les abondances absolues. La méthode d'inférence "mb" consiste à ajuster des régressions pénalisées en utilisant tour à tour chaque espèce comme réponse et toutes les autres comme prédicteurs. De cette manière, le réseau est inféré à partir de 80% des échantillons n fois. Le modèle réalise ces n itérations (inférence du réseau à partir d'un sous-échantillon) sur une plage de valeurs de λ le paramètre qui contrôle la puissance de la régularisation. SpiecEasi génère donc plusieurs tables d’adjacence à partir du même jeux de données initiale, mais de sous-échantillons différents, pour chaque valeur de lambda. Par le biais de l'outil StARS, il sélectionne λ tel que la stabilité globale des arêtes à travers les itérations soit maximisée.La dernière étape consiste à générer le réseau avec cette fois-ci 100% des données et le λ optimisé.
![alt text](https://github.com/paoluxe/How-to-use-SpiecEasi-/blob/ef1049e3acef34d2c039fd7aee9abb8fbc162609/Pictures%20README/Principe%20SpiecEasi.PNG)


## **I) Telechargement des packages**  
SpiecEasi  
```{r eval=FALSE, include=TRUE}
install.packages("devtools")
library(devtools)
install_github("zdk123/SpiecEasi")
```
phyloseq  

```{r eval=FALSE, include=TRUE}
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

```


## **II) Chargement des données et création de l'objet phyloseq**  

SpiecEasi est utilisable avec une matrice d'abondances d'OTUs contenant les échantillons en lignes et les taxons en colonne.
Noms des lignes = Noms des échantillons / noms des colonnes = noms des taxons.

On peut également utiliser un objet phyloseq.

Pour le chargement des données et la création de l'objet phyloseq, je me suis inspirée du lien suivant : https://vaulot.github.io/tutorials/Phyloseq_tutorial.html
Forme et contenu des matrices d'entrée : cf lien ci-dessus
#IMAGE DE AL FORME DES CSV

Importer les données csv

otu_mat = matrice d'abondance des taxons à travers de multiples échantillons
```{r eval=FALSE, include=TRUE}
otu_mat<- read.csv("otu_mat.csv", sep=";", header = T) ; View(otu_mat)
```
tax_mat = matrice de taxinomie
```{r eval=FALSE, include=TRUE}
tax_mat<- read.csv("tax_mat.csv", sep=";", header = T) ; View(tax_mat)
```
samples_df = matrice des variables
```{r eval=FALSE, include=TRUE}
samples_df <- read.csv("samples_df.csv", sep=";", header = T) ; View(samples_df)
```
Ajouter des noms de lignes aux matrices
```{r eval=FALSE, include=TRUE}
library(dplyr)
```
Les noms des lignes deviennent les otu
```{r eval=FALSE, include=TRUE}
otu_mat <- otu_mat %>%
  tibble::column_to_rownames("otu") ## remplacer "otu" par le nom de la colonne contenant les otu/taxons/sp
```
Les noms des lignes deviennent les otu
```{r eval=FALSE, include=TRUE}
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("otu")
```
Les noms des lignes deviennent les échantillons
```{r eval=FALSE, include=TRUE}
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") ## remplacer "sample" par le nom de la colonne contenant les échantillons
```
Transformation des data.frame en matrices 
```{r eval=FALSE, include=TRUE}
otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)
```
Création de l'objet phyloseq
```{r eval=FALSE, include=TRUE}
library(phyloseq)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)

obj_physeq <- phyloseq(OTU, TAX, samples)
```


## **III ) Traitement du jeu de données : selection de groupes, filtre de prevalence/abondance** 
Si tu veux faire plusieurs modalités de réseaux : créer des groupes d'échantillons à analyser
```{r eval=FALSE, include=TRUE}
obj_physeq1 <- subset_samples(obj_physeq, nomvariable =="modalite1") ## remplacer nomvariable et modalite
obj_physeq2 <- subset_samples(obj_physeq, nomvariable =="modalite2") ## remplacer nomvariable et modalite
```

De l'importance des filtres de prévalence
Röttjers et Faust (2018) ont recommandé l'application d'un filtre de prévalence sur les matrices d'abondance, afin d'éviter l'inférence d'associations biaisées, due au nombre élevé de zeros dans les tables d'abondances microbienne et à un effet préférence de niche.

En effet, dans un environnement hétérogène (condition biotiques et abiotiques changeantes entre les échantillons) 2 espèces ayant des optima de croissance dans les mêmes niches co-occureront. L'effet inverse entraînera une exclusion mutuelle.

Les tables d’abondances microbiennes sont riches en zeros, mais on ne connait pas leur signification, est-ce du à une véritable absence ou un sous-échantillonnage?
Les corrélations/régressions calculées sur de nombreux zéros correspondants seront fortement significatives, même si les taxons concernés peuvent varier de façon aléatoire en dessous de la limite de détection. 

si tu veux appliquer un filtre d'abondance ou de prévalence ou les deux arbitrairement
cette ligne signifie : conserver les taxons présent au moins une fois dans 20% des échantillon
ici tu peux jouer sur l'abondance via "sum(x > 0)" et sur la prévalence via  "(0.2*length(x))".
```{r eval=FALSE, include=TRUE}
obj_physeq1_filtre = filter_taxa(obj_physeq1, function(x) sum(x > 0) > (0.2*length(x)), TRUE)
```
tu peux également filtrer les données de façon moins arbitraire via l'outil multicola.
Lien vers l'article : https://academic.oup.com/nar/article/38/15/e155/2409766?login=true
doi:10.1093/nar/gkq545
scripts disponibles ici : https://www.mpi-bremen.de/en/Softwares.html#section1550

## **IV) Réseaux de co-occurrence**  

Obtention du réseau : valeurs des paramètres conseillées par Kurtz
Pour te donner une idée du temps requis, une matrice de 100 échantillons et 400 OTU --> 20 minutes
une matrice de 100 échantillons et 1400 OTU --> 40 minutes
```{r eval=FALSE, include=TRUE}
library(SpiecEasi)
se = spiec.easi(obj_physeq1_filtre,method = 'mb',
                     lambda.min.ratio = 1e-2,
                     nlambda = 20,
                     icov.select.params = list(rep.num = 50))
```
Matrice adjacente "qualitative", une interaction entre deux noeud prend la valeur 1,
quelque soit sa nature, positive ou négative.
```{r eval=FALSE, include=TRUE}
refit_matrix = as.matrix(getRefit(se))
colnames(refit_matrix) <- rownames(otu_table(obj_physeq))
rownames(refit_matrix) <- rownames(otu_table(obj_physeq))
View(refit_matrix)
```
Matrice adjacente "quantitative", une interaction entre deux noeud est pondérée par un coefficient
qui indique la force de l'interaction. Les valeurs vont de -1 à 1
```{r eval=FALSE, include=TRUE}
optbeta_matrix = as.matrix(getOptBeta(se.data))
colnames(optbeta_matrix) <- rownames(otu_table(obj_physeq))
rownames(optbeta_matrix) <- rownames(otu_table(obj_physeq))
```
Matrice de stabilité des arêtes : une interaction entre deux noeud est pondérée par une valeur de stabilité
de l'arête. La stabilité globale des arêtes est calculée à travers les itérations pour le lambda optimale.
```{r eval=FALSE, include=TRUE}
optmerge_matrix = as.matrix(getOptMerge(se.data))
colnames(optmerge_matrix) <- rownames(otu_table(obj_physeq))
rownames(optmerge_matrix) <- rownames(otu_table(obj_physeq))
```
Tu es sensée obtenir les mêmes arêtes avec refit_matrix et optbeta_matrix, ce sont juste
leurs valeurs qui changent, ces arêtes sont calculées à la fin de la procédure avec 100% des échantillons et le 
lambda optimale.

optmerge_matrix te donne la stabilité des arêtes au travers des N répétition (ici 50) pour le lambda optimale
tu n'aura donc pas le meme nombre d'aretes sur les deux précédentes matrices. 
Tu peux tres bien choisir de ne regarder ta matrice que pour une stabilité superieure à un seuil.

Pour une explication (succinte) des fonctions : 
```{r eval=FALSE, include=TRUE}
?getRefit
```
tu verra que les trois que je t'ai proposé ne marchent pas forcément avec un graphe
inféré par la méthode "glasso".

## **V) Calcule de métriques à l'échelle du réseau** 
à partir de liens non quantitatifs
```{r eval=FALSE, include=TRUE}
library(igraph)
ecount(refit_matrix)
average_path_length(refit_matrix)
```

## **VI) Calcule de métriques à l'échelle du noeud** 
à partir de liens non quantitatifs
```{r eval=FALSE, include=TRUE}
library(igraph)
betweenness(refit_matrix, normalized=TRUE)
closeness(refit_matrix, normalized=TRUE)
degree(refit_matrix, normalized=TRUE)
transitivity(refit_matrix, type="local")

```

