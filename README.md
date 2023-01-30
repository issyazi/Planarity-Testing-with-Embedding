### МИНИСТЕРСТВО НАУКИ И ВЫСШЕГО ОБРАЗОВАНИЯ РОССИЙСКОЙ ФЕДЕРАЦИИ <br>
### Федеральное государственное автономное образовательное учреждение высшего образования <br>
## Дальневосточный федеральный университет <br>

### Институт математики и компьютерных технологий
### Департамент информационных и компьютерных систем
### Отчёт о практическом задании по предмету АИСД

# Укладка планарного графа
### Просина Анна Алексеевна, гр. Б9121-09.03.03пикд

<h3 align="center">г. Владивосток, 2022 г.  </h3>
<hr>

# Содержание

- [*Содержание*](#содержание)

- [*Глоссарий*](#глоссарий)

- [*1 Введение*](#введение)

    - [*1.1 Неформальная постановка задачи*](#неформальная-постановка-задачи)

    - [*1.2 Применение*](#применение)

    - [*1.3 Обзор существующих методов*](#обзор-существующих-методов)

- [*2 Гамма-алгоритм*](#гамма-алгоритм)

    - [*2.1 Описание алгоритма*](#описание-алгоритма)

        - [*2.1.1 Входные данные*](#входные-данные)

        - [*2.1.2 Инициализация алгоритма*](#инициализация-алгоритма)

        - [*2.1.3 Общий шаг*](#общий-шаг)

        - [*2.1.4 Завершение работы*](#завершение-работы)


- [*Список литературы*](#список-литературы)


# Глоссарий
- Гомеоморфизм - два графа гомеоморфны (или тождественны с точностью до вершин степени 2), если они оба могут быть получены из одного и того же графа «включением» в его ребра новых вершин степени 2. Если ребро графа изображено в виде линии, то можно на ней поставить точку и считать ее новой вершиной степени 2. 

- Граф – это топологичекая модель, которая состоит из множества вершин и множества соединяющих их рёбер. При этом значение имеет только сам факт, какая вершина с какой соединена.

- Инцидентность (смежность) – отношение между двумя вершинами, в котором существует ребро их соединяющее.

- Мост - ребро, после удаления которого, граф распадается на две компоненты связности. 

- Планарность - свойство графа, которое возникает, если существует укладка взятого графа на плоскости.

- Плоский граф - граф, изображенный на плоскости так, что никакие два ребра не имеют общих точек, кроме инцидентной им обоим вершины.

- Подграф - это часть графа, в которой мы берем некоторые его вершины и ребра. Другими словами, граф H является подграфом графа G, если вершины и ребра H являются подмножеством вершин и ребер G.

- Связный граф - граф, содержащий ровно одну компоненту связности. Это означает, что между любой парой вершин этого графа существует как минимум один путь.

- Сегмент - компонент связности графа с его подграфом.

- Укладка графа - такое его геометрическое изображение, при котором ребра графа пересекаются только в вершинах.


# Введение
Планарная укладка графа (planarity testing with embedding) - это алгоритм проверки графа на планарность и его укладки на плоскость.
Результат работы алгоритма - граф, уложенный без пересечения его рёбер.

## Неформальная постановка задачи
Цель данной работы:
- Изучить по литературным источникам данную задачу и описать ее в форме научного доклада
- Реализовать один из алгоритмов для решения задачи с использованием языка программирования С++
- Исследовать работу алгоритма и описать её

## Применение 
Примером применения данного алгоритма может послужить проблема изготовления электронных микросхем. Электрические цепи печатным способом наносятся на плату из изолирующего материала. Так как наносимые цепи не изолированы, то они не должны пересекаться. В контексте данной проблемы, важно решить, как расположить контакты на схеме, чтобы можно было без пересечений нанести цепи на плату. 

## Обзор существующих методов
1. Теорема Понтрягина-Куратовского - доказывает, что граф планарен тогда и только тогда, когда он не содержит подграфов, гомеоморфных K5 или K3,3. Так как этот критерий очень трудно проверить на практике, данная теорема представляет лишь теоретический интерес.
2. Гамма алгоритм – алгоритм, основанный на теореме Понтрягина-Куратовского, с помощью которого удобнее всего проверить граф на планарность и уложить его на плоскости. 
3. The Hopcroft-Tarjan Planarity Algorithm – алгоритм который проверяет граф на планарность и в зависимости от результата проверки также укладывает его на плоскости

# Гамма-алгоритм
Одним из варинатов, чтобы проверить планарность графа и произвести его плоскую укладку является гамма-алгоритм. Он был подробно описан Ириневым Антоном и Кашириным Виктором в работе под названием “Алгоритм плоской укладки графов”. 

## Описание-алгоритма
### Входные данные
На вход алгоритму подаются графы со следующими свойствами:
1.	Граф связный.
> Если нарушено свойство, то граф нужно укладывать отдельно по компонентам связности. 

2.	Граф содержит хотя бы один цикл.
> Если нарушено свойство , то граф — дерево и нарисовать его плоскую укладку тривиально.

3.	Граф не имеет мостов.
> Случай нарушения данного свойства рассмотрим более подробно. Если в графе есть мосты, то их нужно разрезать, провести отдельно плоскую укладку каждой компоненты связности, а затем соединить их мостами. Здесь может возникнуть трудность: в процессе укладки концевые вершины моста могут оказаться внутри плоского графа. Нарисуем одну компоненту связности, и будем присоединять к ней другие последовательно. Каждую новую компоненту связности будем рисовать в той грани, в которой лежит концевая вершина соответствующего моста. Так как граф связности мостами компонент связности является деревом, мы сумеем получить плоскую укладку.

**Замечание.** Если концевая вершина, принадлежащая новой (только что нарисованной) компоненте связности, также оказалась внутри, необходимо выполнить повторную укладку этой компоненты так, чтобы ребро, содержащее концевую вершину, принадлежало внешней грани.

### Инициализация алгоритма
<p align="center">
<image
  src="/images/original_graph.svg"
  alt="Исходный граф"
  caption="Исходный граф"
  style="width: 400px;">
  <p align="center">Исходный граф</p>
</p>

выбираем любой простой цикл в G, пусть это будет {1, 2, 3, 4, 5, 6}, укладываем его на плоскости, и получаем две грани: Г1 — внешнюю и Г2 — внутреннюю (см. рис. 5).

Уже уложенную часть исходного графа будем обозначать как G′, тогда после первого шага G′ представляет собой цикл {1, 2, 3, 4, 5, 6}.

На каждом шаге будем строить множество сегментов. Каждый сегмент S относительно уже построенного графа G′ представляет собой одно из двух:

- ребро, оба конца которого принадлежат G′, но само оно не принадлежит G′;
- связную компоненту графа G – G′, дополненную всеми ребрами графа G, один из концов которых принадлежит связной компоненте, а второй из графа G′.

Вершины, которые одновременно принадлежат G′ и какому-то сегменту, назовем контактными вершинами. Для нашего примера сегменты изображены на рис. 6. Контактные вершины обведены в квадрат.

Если бы в каком-нибудь сегменте не было ни одной контактной вершины, то граф до построения множества сегментов был бы несвязный; если бы была только одна контактная вершина, то граф имел бы мост. Эти возможности заранее исключены, так что каждый сегмент имеет не менее двух контактных вершин. Поэтому в каждом сегменте имеется цепь между любой парой таких вершин.

Если все контактные вершины сегмента S имеют номера вершин какой-то грани Г, то мы будем говорить, что грань Г вмещает этот сегмент, в этом случае будем использовать следующее обозначение: S Г. Однако, может быть так, что не одна грань вмещает в себя сегмент S, а несколько. Множество таких граней обозначим Г(S), а их число |Г(S)|.

### Общий шаг 
Выделяются все сегменты Si и определяются числа |Г(Si)|. Если хоть одно из них равно 0, то граф не планарен, конец. Иначе, выбираем сегмент, для которого число |Г(S)| минимально, или любой из них, если таких сегментов несколько. В этом сегменте найдем произвольную цепь между двумя контактными вершинами и уложим ее в любую из граней множества Г(S). При этом данная грань разобьется на две. Уже уложенная часть графа G’ после укладки цепи увеличится, а сегмент, из которого вынута цепь, исчезнет или развалится на меньшие с новыми контактными вершинами, ведущими к вершинам G′.

В результате повторения общего шага будет либо получена плоская укладка, когда множество сегментов станет пустым, либо будет получено, что граф G не является планарным.

Вернемся к нашему примеру. Пока для любого i: Si {Г1, Г2}, |Г(Si)| = 2. Поэтому возьмем первый по номеру сегмент Si и в нем цепь {1, 4}; вставим эту цепь в грань Г2. После укладки цепи G’ увеличится и произойдут изменения в структуре сегментов (см. рис. 7 a, b).

Определим, какие грани вмещают новые сегменты. Теперь сегменты S1 и S3 можно уложить только в одну грань Г1, в то время как сегменты S2 и S4 можно уложить в две грани (для S2 это грани Г1 и Г2, для S4 - Г1 и Г3). Поэтому берем S1. Возьмем в нем цепь {2, 5} и уложим ее в Г1. Получим увеличенный граф G′ и уменьшенную систему сегментов (см. рис. 8 a, b).


### Завершение работы

Продолжая таким образом, в итоге получим плоскую укладку исходного графа G.

# Список литературы
1. [Graph Planarity Testing with Hierarchical Embedding Constraints Giuseppe Liottaa, Ignaz Rutterb, Alessandra Tappinia.](https://arxiv.org/pdf/1904.12596.pdf)
2. [Иринёв Антон, Каширин Виктор Алгоритм плоской укладки графов.](https://studfile.net/preview/1869599/)
3. [Kurt Mehlhorn and Petra Mutzel., On the embedding phase of the Hopcroft and Tarjan planarity testing algorithm., Algorithmica, 1996](https://domino.mpi-inf.mpg.de/internet/reports.nsf/efc044f1568a0058c125642e0064c817/02b4941bb1079240c12560b700590d27/$FILE/MPI-I-94-117.pdf)
4. [Статья на habr.com. Теория Графов.](https://habr.com/ru/post/565998/)
5. https://neerc.ifmo.ru/wiki/index.php?title=Гамма-алгоритм
6. https://neerc.ifmo.ru/wiki/index.php?title=Укладка_графа_на_плоскости – укладка графа на плоскости
7. https://neerc.ifmo.ru/wiki/index.php?title=Теорема_Понтрягина-Куратовского – теорема, которая доказывается с помощью гамма-алгоритма
8. https://portal.tpu.ru/SHARED/t/TRACEY/Courses/Graph_Theory/Tab1/graph_lec_09.pdf - алгоритм укладки планарного графа на плоскости
9. https://github.com/tehnik819/Gamma-Algorithm-Java/blob/master/src/Graph.java - код
10. https://www.youtube.com/watch?v=MS98VuAG9Yo – как нарисовать плоский граф через гамма-алгоритм
11. https://users.math-cs.spbu.ru/~okhotin/teaching/tcs1_2016/okhotin_tcs1_2016_l5.pdf - лекция по плоским графам
12. https://dic.academic.ru/dic.nsf/ruwiki/628470 - гамма-алгоритм
13. https://textarchive.ru/c-2330757-p3.html - 
14.	https://theslide.ru/uncategorized/osnovnye-ponyatiya-teorii-grafov - теория графов
15.	https://studopedia.ru/2_84961_ploskie-i-planarnie-grafi-ploskie-karti-teorema-eylera.html - плоские и планарные графы
16.	https://progr-system.ru/wp-content/uploads/Math/МАПКС-08-ПлоскиеПланарныеГрафыРаскраска.pdf - плоские и планарные графы
17.	https://prezi.com/p/ezmde7t_sf6h/presentation/ - алгоритмы проверки планарности графов
18.	https://www.researchgate.net/publication/228898512_The_Gamma_Algorithm_in_Convex_Cone_Analysis_of_Hyperspectral_Images - гамма алгоритм
19.	https://www.youtube.com/watch?v=das1x_ntedU - Укладка планарных графов, алгоритм D. Eppstein
20.	 https://studfile.net/preview/1869599/ - гамма алгоритм для плоской укладки графа
21.	https://cs.brown.edu/people/rtamassi/gdhandbook/chapters/planarity.pdf
22.	https://www.researchgate.net/publication/228993896_The_Hopcroft-Tarjan_Planarity_Algorithm - укладка планарного графа
23.	https://en.wikipedia.org/wiki/Planarity_testing 
24.	https://towardsdatascience.com/graph-planarity-and-path-addition-method-of-hopcroft-tarjan-for-planarity-testing-c56d2df2b0b3 - hopcroft-tarjan planarity testing 
25.	https://github.com/shawnwanderson/Hopcroft-Tarjan-Planarity-Testing - code c++ of Hopcroft tarjan testing of planarity
26.	 https://jgaa.info/accepted/2008/GutwengerKleinMutzel2008.12.1.pdf
27.	https://link.springer.com/content/pdf/10.1007/978-3-540-77537-9_9.pdf 
