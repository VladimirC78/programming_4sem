# Лабораторная работа №3

В этой работе решается уравнение Пуассона с измененными краевыми условиями и областью задания. 

## Файлы проекта

domain_creator.cpp создает область задания уравнения (пятиконечная звезда) и сетку на ней, сохраняет ее в формате .msh.

main.cpp решает уравнение Пуассона с заданными граничными условиями и выводит результат в папку snapshots.

## Сборка и запуск

Для создания заголовочного файла Poisson.h введите в терминал в папке Lab03:

``` bash
ffc -l dolfin Poisson.ufl
```

Для запуска выполнить

``` bash
cmake -B build
cmake --build build
build/demo_poisson
build/domain_creator
```

## Получение сетки в формате .xml

Сетку в формат .xml можно преобразовать с помощью программы meshio, для этого нужно выполнить:

``` bash
meshio convert star.msh domain.xml
```

И в файле .xml поставить dim="2" вместо "3"!

## Результат

![image](https://github.com/VladimirC78/programming_4sem/blob/master/Lab03/star_with_heart.png)
