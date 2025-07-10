#!/bin/bash

# Задайте путь к директории
DIRECTORY="/mnt/ligandpro/data/garifullin/ligand_protein_system_preparation/train"  # Замените на нужный путь

# Задайте имя выходного файла
OUTPUT_FILE="MOAD_train"

# Создать или очистить выходной файл
> "$OUTPUT_FILE"

# Перебрать все директории и записать их названия в файл
for dir in "$DIRECTORY"/*/; do
    # Проверяем, является ли это директорией
    if [ -d "$dir" ]; then
        # Получаем только имя папки и добавляем его в выходной файл
        echo "$(basename "$dir")" >> "$OUTPUT_FILE"
    fi
done

echo "Список папок сохранён в $OUTPUT_FILE"
