# ==========================================================================
# Packages used
using PyPlot

# ==========================================================================
# Parameters for the resolution algorithm
srand(1)
seedingStrategy = false # active ou pas la strategie des solutions seedings
nbInd = 50              # nombre d'individus !! doit etre un nombre pair
maxGeneration = 100     # nombre de generation de l'algorithme genetique

# ==========================================================================
# Type of a instance for the scheduling problem presented in the paper Morita et al., 2001 (FCDS)
type t_instance
  n ::Int64      # number of tasks
  p              # processing time
  d              # due date
  w              # weight
end

# ==========================================================================
# Type of a solution for the scheduling problem presented in the paper Morita et al., 2001 (FCDS)
type t_solution
  x               # vector of index of jobs in the solution
  z1    ::Int64   # performance objective 1
  z2    ::Int64   # performance objective 2
end

# ==========================================================================
# Create a didactic instance
function createDidacticInstance()
  tasks = 4
  data = t_instance(tasks, zeros(Int64,tasks), zeros(Int64,tasks), zeros(Int64,tasks))
  data.p = [3, 4, 5, 6]
  data.d = [20, 16, 11, 5]
  data.w = [1, 1, 1, 1]
  return data
end
#data = createDidacticInstance()

# ==========================================================================
# Create a parameterized random instance
function createRandomInstance(tasks, maxp, maxd, maxw)
  data = t_instance(tasks, zeros(Int64,tasks), zeros(Int64,tasks), zeros(Int64,tasks))
  data.p = rand(1:maxp, tasks)
  data.d = rand(1:maxd, tasks)
  data.w = rand(1:maxw, tasks)
  return data
end

# ----------------------------------
# minimize the flowtime on a single machine (Smith's rule)
function minFlowtime(data::t_instance)
  # minimise f1 : flowtime (SPT-rule)
  SPT = sortperm(data.p)
  somCi = 0
  somPi = 0
  c = zeros(Int64,data.n)
  for i=1:data.n
    c[SPT[i]] = somPi + data.p[SPT[i]]
    somPi = somPi + data.p[SPT[i]]
    somCi = somCi + c[SPT[i]]
  end
  println("minimum f1 (min flowtime) = ", somCi)
  return SPT, somCi
end

# ----------------------------------
# minimize the maximum tardiness on a single machine
function minMaxTardiness(data::t_instance)
  # minimise f2 : maximum tardiness (EDD-rule)
  EDD = sortperm(data.d)
  maxT = 0
  somPi = 0
  c = zeros(Int64,data.n)
  for i=1:data.n
    c[EDD[i]] = somPi + data.p[EDD[i]]
    somPi = somPi + data.p[EDD[i]]
    maxT = max(maxT, c[EDD[i]]-data.d[EDD[i]])
  end
  println("minimum f2 (max tardiness) = ", maxT)
  return EDD, maxT
end

# ==========================================================================
# Bi-objective evaluation of a solution
function evaluateSolution(data::t_instance, sol::t_solution)
  # a solution here is a feasible scheduling of tasks
  c = 0; sol.z1 = 0; sol.z2 = 0
  for i=1:data.n
    # completion time de la tache i
    c = c + data.p[sol.x[i]]
    # calcule la somme des c_i (total flow time))
    sol.z1 = sol.z1 + c
    # calcule le retard maximum (maximum tardiness)
    sol.z2 = max(sol.z2, c - data.d[sol.x[i]] )
  end
end

# ==========================================================================
# Generate and evaluate an initial randomized population
function generateInitialPopulation(data::t_instance, nSol::Int64)
  # pop est un tableau de 1..nSol de type solution
  pop = Array{t_solution}(nSol)
  for i=1:nSol
    # sequence = permutation aleatoire de n entiers
    sol = t_solution(randperm(data.n), 0, 0)
    # evaluation de la sequence sur les 2 objectifs
    evaluateSolution(data, sol)
    # sauvegarde des solutions evaluees dans la population
    pop[i] = sol
  end
  return pop
end

# ==========================================================================
# Generate and evaluate an initial randomized population
function generateInitialPopulationPeek(data::t_instance, nSol::Int64, z1opt::t_solution, z2opt::t_solution)
  # pop est un tableau de 1..nSol de type solution
  pop = Array{t_solution}(nSol)
  pop[1] = doMutation(data, z1opt)
  for i=2:floor(Int64,nSol/2)
    pop[i] = doMutation(data, pop[i-1])
  end
  pop[floor(Int64,nSol/2)+1] = doMutation(data, z2opt)
  for i=floor(Int64,nSol/2)+2:nSol
    pop[i] = doMutation(data, pop[i-1])
  end

  return pop
end

# ----------------------------------
# Get the maximal value on both objectives for a given population
function getSupRangesObjectives(pop)
  # range max de l'abscisse et l'ordonnee du nuage
  rangeMaxX = 0
  rangeMaxY = 0
  # extrait les performances des points du nuage
  for i=1:length(pop)
    individu = pop[i]
    rangeMaxX = max(rangeMaxX, individu.z1)
    rangeMaxY = max(rangeMaxY, individu.z2)
  end
  # rend les axes orthonormes
  rangeXY = max(rangeMaxX, rangeMaxY)
  return rangeXY
end

# ----------------------------------
# Plot the performances (z1,z2) of a given population
function plotPerformancesPopulation(pop)
  # vecteurs abscisse et ordonnee du nuage
  x = zeros(Int64,length(pop))
  y = zeros(Int64,length(pop))
  # extrait les performances des points du nuage
  for i=1:length(pop)
    individu = pop[i]
    x[i] = individu.z1
    y[i] = individu.z2
  end
  # trace les points
  plot(x[:],y[:], color="black",linestyle="",marker="+")
end

# ----------------------------------
# Apply on two individuals the crossover operator defined by Muhlenbein et al. 1988
function doCrossover1988(data::t_instance, donator::t_solution, receiver::t_solution)
# Operateur de crossover de Muhlenbein et al. 1988 (FCDS 2001)
  n = length(donator.x) # nbre de taches dans une solution

  # extrait une sequence contigue circulaire de taches du donator
  longueurSequence = floor(Int64, n * rand(100:500)/1000) + 1
  iDeb = rand( 1:(n-longueurSequence+1) )
  childA = donator.x[iDeb:iDeb+longueurSequence-1]

  # extrait du receiver toutes les taches non presentes dans la sequence extraite du donator
  lChildB = n - length(childA)
  j=1
  childB = zeros(Int64,lChildB)
  for i=1:n
    if !(receiver.x[i] in childA)
      childB[j]= receiver.x[i]
      j=j+1
    end
  end
  # construit l'enfant ainsi cree
  child = vcat(childA,childB)

  # evaluation de l'enfant au regard des 2 objectifs
  enfant = t_solution(child, 0, 0)
  evaluateSolution(data, enfant)

  return enfant
end

# ----------------------------------
# Apply on two individuals a crossover operator (modified version of by Muhlenbein et al. 1988)
function doCrossover1988bis(data::t_instance, donator::t_solution, receiver::t_solution)
# Operateur de crossover de Muhlenbein et al. 1988 modifie (extraction criculaire du donateur)
  n = length(donator.x) # nbre de taches dans une solution

  # extrait une sequence contigue circulaire de taches du donator
  iDeb = rand(1:n)
  iFin = iDeb + floor(Int64, n * rand(100:500)/1000)
  lChildA = iFin-iDeb+1
  j=1
  childA = zeros(Int64,lChildA)
  for i=iDeb:iFin
    #println(donator.x[(i%n)+1])
    childA[j]= donator.x[(i%n)+1]
    j=j+1
  end
  # extrait du receiver toutes les taches non presentes dans la sequence extraite du donator
  lChildB = n - lChildA
  j=1
  childB = zeros(Int64,lChildB)
  for i=1:n
    if !(receiver.x[i] in childA)
      childB[j]= receiver.x[i]
      j=j+1
    end
  end
  # construit l'enfant ainsi cree
  child = vcat(childA,childB)

  # evaluation de l'enfant au regard des 2 objectifs
  enfant = t_solution(child, 0, 0)
  evaluateSolution(data, enfant)

  return enfant
end

# ----------------------------------
# Apply on one individual the mutation operator defined by a swap between tasks
function doMutation(data::t_instance, individu::t_solution)
# Operateur de mutation (FCDS 2001)

  # clone l'individu recu comme individu mutant
  mutant = deepcopy(individu)
  mutant.z1 = 0; mutant.z2 = 0

  # extrait deux taches et echange de celles-ci dans la solution
  n = length(mutant.x)                                    # nbre de taches dans une solution
  iTacheA = rand( 1:n );       iTacheB =  rand( 1:n )      # indice des taches a echanger
  tacheA = mutant.x[iTacheA];  tacheB = mutant.x[iTacheB] # taches a echanger
  mutant.x[iTacheB] = tacheA;  mutant.x[iTacheA] = tacheB # echange

  # evaluation du mutant au regard des 2 objectifs
  evaluateSolution(data, mutant)

  return mutant
end

# ==========================================================================
# Compute the fitness of all individuals of a population according the Goldberg definition for 3 cases
function computeGoldbergRank3Cases(pop, cas)
# L'algorithme fonctionne par strates; il prend comme hypothese que tous les points sont de rang 1.
# Si ce n'est pas le cas, il corrige le rang et recommance avec les rangs de plus grandes valeurs.

# Deux fonctions sont a minimiser avec 3 cas : cas 1 = S.O.; cas 2 = N.O.; cas 3 = S.E.

  # vecteur qui memorise le rang des points de la population
  taillePop = length(pop)
  rang = ones(Int64, taillePop)
  rangCourant = 0
  while true # repeat...

    # les rangs de niveaux inferieurs au niveau courant sont figes.
    # Tous les autres sont alignes sur le rang courant.
    calculRangIncomplet = false
    rangCourant =  rangCourant +1
    for i=1:taillePop
      if (rang[i] >  rangCourant)
        rang[i] = rangCourant
      end
    end
    for i=1:taillePop
      for j=1:taillePop
        if cas == 1 #---------------------------------------
          # minZ1 minZ2
          if ( (rang[i]>= rangCourant) && (rang[j]>= rangCourant) &&  # individus de rangs libres seulement
               ( (pop[i].z1 > pop[j].z1 && pop[i].z2 > pop[j].z2) ||      # dominance stricte de j sur i
                 (pop[i].z1 == pop[j].z1 && pop[i].z2 > pop[j].z2) ||     # dominance faible de j sur i
                 (pop[i].z1 > pop[j].z1 && pop[i].z2 == pop[j].z2)        # dominance faible de j sur i
               )
             )
            rang[i]=rang[i]+1
            calculRangIncomplet = true
          end
        elseif cas == 2 #---------------------------------------
          # minZ1 maxZ2
          if ( (rang[i]>= rangCourant) && (rang[j]>= rangCourant) &&
               ( (pop[i].z1 > pop[j].z1 && pop[i].z2 < pop[j].z2) ||
                 (pop[i].z1 == pop[j].z1 && pop[i].z2 < pop[j].z2) ||
                 (pop[i].z1 > pop[j].z1 && pop[i].z2 == pop[j].z2)
               )
             )
            rang[i]=rang[i]+1
            calculRangIncomplet = true
          end
        else # cas==3 ---------------------------------------
          # maxZ1 minZ2
          if ( (rang[i]>= rangCourant) && (rang[j]>= rangCourant) &&
               ( (pop[i].z1 < pop[j].z1 && pop[i].z2 > pop[j].z2) ||
                 (pop[i].z1 == pop[j].z1 && pop[i].z2 > pop[j].z2) ||
                 (pop[i].z1 < pop[j].z1 && pop[i].z2 == pop[j].z2)
               )
             )
            rang[i]=rang[i]+1
            calculRangIncomplet = true
          end
        end
      end #for j
    end #for i
    calculRangIncomplet == false && break # ...until
  end
  return rang
end

# ==========================================================================
# Compute the fitness of all individuals of a population according when both objectives have to be minimized
function computeFitnessPopulation(pop)

  # Calcule le rang d'une population selon la definiton de Goldberg dans 3 secteurs du plan (SO; NO; SE)
  cas=1; rangSudOuest  = computeGoldbergRank3Cases(pop, cas) # SO
  cas=2; rangNordOuest = computeGoldbergRank3Cases(pop, cas) # NO
  cas=3; rangSudEst    = computeGoldbergRank3Cases(pop, cas) # SE

  # calcule un point median de reference deduit des points maximaux sur les 2 objectifs
  # passage par vZ par pas trouve comment faire directement sur pop
  #vZ = zeros(Int64,length(pop),2)
  #for i=1:length(pop) vZ[i,1] = pop[i].z1 ; vZ[i,2] = pop[i].z2 end
  #iMaxZ1 = indmax(vZ[:,1]) ; iMaxZ2 = indmax(vZ[:,2])

  #z1Limit = floor(Int64, (vZ[indmax(vZ[:,1]),1] + vZ[indmax(vZ[:,2]),1])/2)
  #z2Limit = floor(Int64, (vZ[indmax(vZ[:,1]),2] + vZ[indmax(vZ[:,2]),2])/2)

  zIdeal = zeros(Int64,2)
  zAntiIdeal = zeros(Int64,2)
  zNadir = zeros(Int64,2)
  vZ = zeros(Int64,length(pop),2)
  for i=1:length(pop) vZ[i,1] = pop[i].z1; vZ[i,2] = pop[i].z2 end
  iMaxZ1 = indmax(vZ[:,1]) ; iMaxZ2 = indmax(vZ[:,2])
  iMinZ1 = indmin(vZ[:,1]) ; iMinZ2 = indmin(vZ[:,2])

  zIdeal[1]=vZ[iMinZ1,1] ; zIdeal[2]=vZ[iMinZ2,2]
  zAntiIdeal[1]=vZ[iMaxZ1,1] ; zAntiIdeal[2]=vZ[iMaxZ2,2]
  zNadir[1]=vZ[iMinZ2,1] ; zNadir[2]=vZ[iMinZ1,2]
#  z1Limit=zNadir[1] ; z2Limit=zNadir[2]
  z1Limit = floor(Int64, (zIdeal[1] + zAntiIdeal[1])/2 - 0.025* zAntiIdeal[1])
  z2Limit = floor(Int64, (zIdeal[2] + zAntiIdeal[2])/2 - 0.025* zAntiIdeal[2])

  plot(zIdeal[1],zIdeal[2], color="red",linestyle="",marker="+")
  plot(zAntiIdeal[1],zAntiIdeal[2], color="red",linestyle="",marker="+")
  plot(z1Limit,z2Limit, color="red",linestyle="",marker="*")

#  plot(zNadir[1],zNadir[2], color="red",linestyle="",marker="*")

  # neutralise les rangs des individus domines par z_ref pour les secteurs NO et SE
  for i=1:length(pop)
    if (vZ[i,1] >z1Limit) && (vZ[i,2] >z2Limit)
      rangNordOuest[i] = length(pop) # comme valeur infinie
      rangSudEst[i]    = length(pop) # comme valeur infinie
    end
  end

  # neutralise les rangs > limite des individus pour les secteurs NO et SE
  for i=1:length(pop)
    if rangNordOuest[i] > 1 #length(pop)
      rangNordOuest[i] = length(pop) # comme valeur infinie
    end
    if rangSudEst[i] > 1# length(pop)
      rangSudEst[i]    = length(pop) # comme valeur infinie
    end
  end

  # calcule les lignes de rangs
#  rangParents  =   rangSudOuest
  rangParents  =   min(rangSudOuest, 2*rangNordOuest, 2*rangSudEst)

  return rangParents, vZ, iMaxZ1, iMaxZ2, z1Limit, z2Limit
end

# ==========================================================================
# Select two individuals by tournament from a given population
function SelectTwoGoodIndividuals(pop, rangPop)

  p1 = 0; p2 = 0; p3 = 0; p4 = 0  # indice des 4 parents en competition
  i1 = 0; i2 = 0                  # indice des 2 parents retenus

  # selectionne par tournoi un premier parent parmis p1 et p2 avec p1 <> p2
  p1 = rand(1:length(pop))
  while true
    p2 = rand(1:length(pop))
    p1 != p2 && break
  end
  if rangPop[p1] == rangPop[p2]
    i1 = rand() < 0.5 ? p1 : p2  # tirage au hazard en cas d'ex-aequo
  else
    i1 = rangPop[p1] < rangPop[p2] ? p1 : p2 # parent de meilleur fitness retenu
  end

  # selectionne par tournoi un second parent parmis p3 et p4 avec p3 <> p4
  p3 = rand(1:length(pop))
  while true
    p4 = rand(1:length(pop))
    p3 != p4 && break
  end
  if rangPop[p3] == rangPop[p4]
    i2 = rand() < 0.5 ? p3 : p4  # tirage au hazard en cas d'ex-aequo
  else
    i2 = rangPop[p3] < rangPop[p4] ? p3 : p4  # parent de meilleur fitness retenu
  end

  # retourne les indices des 2 parents retenus a l'issue du tournoi
  return i1, i2
end

# ==================================
# ==================================
# Met en place une instance numerique
tasks = 100  # nombre de taches dans le probleme
maxp  = 15  # plus grand processing time
maxd  = 100  # plus grande due date
maxw  = 1   # plus grande penalite

# soit 1) instance didactique du papier de 2001
#data = createDidacticInstance()
# soit 2) instance aleatoire parametree
data  = createRandomInstance(tasks, maxp, maxd, maxw)

# ----------------------------------
# Elements pour la representation graphique de l'espace des objectifs
xlabel("z1 (total flow time)")
ylabel("z2 (Maximum tardiness)")

# ==================================
# RESOLUTION EXACTE
# ----------------------------------
# Calcul des solutions optimales pour les 2 objectifs pris separement
optf1= t_solution(zeros(data.n), 0, 0)
optf1.x, optf1.z1 = minFlowtime(data)
evaluateSolution(data, optf1) # (re)evalue la solution sur les 2 objectifs

optf2= t_solution(zeros(data.n), 0, 0)
optf2.x, optf2.z2 = minMaxTardiness(data)
evaluateSolution(data, optf2) # (re)evalue la solution sur les 2 objectifs

plot(optf1.z1,optf1.z2, color="green",linestyle="",marker="D")
plot(optf2.z1,optf2.z2, color="green",linestyle="",marker="D")

# ==================================
# RESOLUTION APPROCHEE
# ----------------------------------
# Population initiale
#pop = generateInitialPopulation(data, nbInd)
pop = generateInitialPopulationPeek(data, nbInd, optf1, optf2)

# strategie de seeding solutions
if seedingStrategy == true
  # insertion de 2 seeding solutions
  pop[1] = deepcopy(optf1)
  pop[2] = deepcopy(optf2)
end

# ----------------------------------
# Elements pour la representation graphique de l'espace des objectifs
rangeXY = getSupRangesObjectives(pop)
#xlim(0,rangeXY+10)
#ylim(0,rangeXY+10)

# ----------------------------------
# Trace la population initiale
plotPerformancesPopulation(pop)

# ----------------------------------
# Calcule le rand des individus de la population
rangPopulation, vZ, iMaxZ1, iMaxZ2, z1Limit, z2Limit = computeFitnessPopulation(pop)

#plot(vZ[iMaxZ1,1],vZ[iMaxZ1,2], color="red",linestyle="",marker="+")
#plot(vZ[iMaxZ2,1],vZ[iMaxZ2,2], color="red",linestyle="",marker="+")
#plot(z1Limit,z2Limit, color="red",linestyle="",marker="*")

# ----------------------------------
# Operations evolutionaires

for generation = 1:maxGeneration
  newPop = []
  for evolutionNaturelle=1:length(pop)/2

    # Selectionne deux bons individus
    i1, i2 = SelectTwoGoodIndividuals(pop, rangPopulation)

    # Crossover entre deux parents et creation d'un enfant
    individuC = doCrossover1988(data, pop[i1],pop[i2])
    push!(newPop,individuC)

    # Mutation de l'enfant
    individuM = doMutation(data, individuC)
    push!(newPop,individuM)

    plot(individuC.z1,individuC.z2, color="blue",linestyle="",marker="o")
    plot(individuM.z1,individuM.z2, color="blue",linestyle="",marker="x")

  end
  #pop = selectNextGeneration(pop, newPop)
  pop = vcat(pop, newPop)
  rangPopulation, vZ, iMaxZ1, iMaxZ2, z1Limit, z2Limit = computeFitnessPopulation(pop)

  # ==========================================================================
  # newgeneration est un tableau de 1..nSol de type solution
  newgeneration = [] #Array{t_solution}(nSol*2)

  els = 1
  rangCourant = 1
  while els <= nbInd # tant que la nouvelle generation n'est pas complete faire
    for i=1:2*nbInd # extrait de la population tous les individus de rang = rangCourant
      if (rangPopulation[i] == rangCourant)
        push!(newgeneration, pop[i])
        rangPopulation[i] =  0
        els = els + 1
      end
    end
    rangCourant = rangCourant +1
  end

  pop = copy(newgeneration[1:nbInd])
  rangPopulation, vZ, iMaxZ1, iMaxZ2, z1Limit, z2Limit = computeFitnessPopulation(pop)
  #return newgeneration[1:nbInd]

  # extrait les performances des points du nuage
  for i=1:nbInd
    plot(pop[i].z1,pop[i].z2, color="orange",linestyle="",marker=".")
  end
  # trace les points
  println("Generation = ", generation)

end

# ==========================================================================
# Extraction des solutions potentiellement non dominees (cad de rang 1 SO)
cas=1; rangSudOuest  = computeGoldbergRank3Cases(pop, cas)
for i=1:nbInd
  if rangSudOuest[i]==1
# if rangPopulation[i]==1
    println(i, ") ", pop[i].z1, "  ",pop[i].z2)
    plot(pop[i].z1,pop[i].z2, color="cyan",linestyle="",marker="D")
  end
end

zIdeal = zeros(Int64,2)
zAntiIdeal = zeros(Int64,2)
zNadir = zeros(Int64,2)
vZ = zeros(Int64,length(pop),2)
for i=1:length(pop) vZ[i,1] = pop[i].z1; vZ[i,2] = pop[i].z2 end
iMaxZ1 = indmax(vZ[:,1]) ; iMaxZ2 = indmax(vZ[:,2])
iMinZ1 = indmin(vZ[:,1]) ; iMinZ2 = indmin(vZ[:,2])

zIdeal[1]=vZ[iMinZ1,1] ; zIdeal[2]=vZ[iMinZ2,2]
zAntiIdeal[1]=vZ[iMaxZ1,1] ; zAntiIdeal[2]=vZ[iMaxZ2,2]
zNadir[1]=vZ[iMinZ2,1] ; zNadir[2]=vZ[iMinZ1,2]
# ----------------------------------
# ----------------------------------
# ----------------------------------
