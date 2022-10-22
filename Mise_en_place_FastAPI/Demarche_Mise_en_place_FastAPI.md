![A cat](C:\Users\cosne\Documents\test_github\projet-aln\Mise_en_place_FastAPI/logo.jpg)



>  >   >   >    > ## Année universitaire 2022-2023
> ## Projet ingénieur 

## Démarche 
## Mise en place FastAPI




Elèves : Cosnelle DJOUMEKOUM & Julie RIDOLFI   
*Diplôme d’ingénieur “Informatique et Statistique” - 5ème année*   


Tuteur école : Frédéric HOOGSTOEL     
*Enseignant universitaire à Polytech’Lille*   

Tuteur entreprise : Damien MARCHAL     
*Ingénieur de recherche CNRS*   






POLYTECH LILLE
Av. Paul Langevin - Cité Scientifique
59650 Villeneuve d’Ascq
(33) 03 28 76 73 60
(33) 03 28 76 73 61


CRISTAL
Av. Henri Poincaré - Cité Scientifique
59650 Villeneuve d’Ascq

&nbsp;

&nbsp;

**FastAPI** est un framework web haute performance, open source, permettant de créer des APIs avec    
Python à partir de la version 3.6          

**Prérequis :**          
Installer  python          

**Démarche mise en place :**         
1.	Se placer dans un terminal et exécuter les commande suivantes :    

    - **pip install fastapi** : permet d’installer FastAPI   
    - **pip install uvicorn** : permet le lancer le serveur API en local   

2.	Créer ensuite son répertoire de travail et y ajouter un fichier .py dans le lequel le code python sera écrit  
Un exemple de fichier  **‘’main.py ‘’** contenant le code pour importer FastAPI et créer une route pour  obtenir nos données au format .json  depuis notre serveur API sera joint à ce document.  

3.	Se placer ensuite dans le répertoire de travail et lancer la commande suivante pour lancer le serveur FastAPI en local :  
**uvicorn main:app --reload**      
       On peut avoir une erreur de reconnaissance de la commande uvicorn par le terminal si on                  effectue l’installation sur windows, dans ce cas utiliser la commande :  
              **python -m uvicorn main:app --reload**  
        ou travailler dans un environnement de travail python :         
              **python -m venv venv** : permet de créer un environnement de travail virtuel nommé venv   
              **venv\Scripts\activate** : permet d’activer l’environnement de travail virtuel   
              **deactivate** :  permet de désactiver l’environnement de travail   
        l’environnement de travail permet également d’éviter de télécharger les modules nécessaires     
        pour un projet directement sur son PC   

Le lien du serveur local sur lequel l’API est lancer est ensuite fourni : http://127.0.0.1:8000   


**A savoir :**   
http://127.0.0.1:8000/docs permet d’avoir de la documentation sur notre API   









