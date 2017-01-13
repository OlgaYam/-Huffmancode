# -Huffmancode
В программе необходимо задать алфавит в виде строки, без пробелов. Далее можно сгенерировать матрицу переходных вероятностей рандомно или считать ее из файла. Затем производится сжатие тремя способами.

Размер текста 100000 символов (параметр size).
При сжатии выхода Марковского источника в предположении, что символы получены из дискретного стационарного источника без памяти, получены следующие результаты:
Теоритическая энтропия = 1.51336
Практическая энтропия = 1.562

При сжатии выход Марковского источника, оценив матрицу переходных вероятностей по выходу источника:
Практическая энтропия = 1.51032

При сжатии выхода Марковского источника, используя истинную матрицу переходных вероятностей:
Практическая энтропия = 1.51032
