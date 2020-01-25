import csv


def read_txt(filename):
    """
    Read file with given filename as txt file and extract rows, strip and clean up, convert and return nested list of
    rows and values

    :param filename:
    :return:
    """
    rows = []

    with open(filename, 'r') as file:
        contents = file.readlines()

        for line in contents:

            line = line[1:]
            line = line[:-2]

            elements = line.split(', ')

            for i in range(len(elements)):
                if elements[i] != 'None':
                    elements[i] = elements[i][1:]
                    elements[i] = elements[i][:-1]

            rows.append(elements)

    return rows


def add_dates_to_list(input_list):
    d = {}

    j = 0
    for year in range(2017, 2020):

        for month in range(1, 13):

            days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

            for day in range(days[month - 1] + 1):
                date = str(year) + '-' + str(month) + '-' + str(day)
                d[date] = input_list[j]
                j += 1

    return d


def farenheit_to_celsius(name):
    """
    Change farenheit to celsius in the csv file

    :param name:
    :return:
    """
    with open(name, 'r', newline='') as csv_file:

        reader = csv.reader(csv_file)

        for row in reader:

            for i in range(len(row)):
                if (i - 1) % 12 == 0 or (i - 2) % 12 == 0:

                    if row[i] != 'None' and row[i] != '--':
                        row[i] = round((float(row[i]) - 32) * 5 / 9, 2)


def str_to_num(lst):
    """
    Convert the entries in the list that are numbers to type float (all except entries i=0,4)

    :param lst:
    :return:
    """
    for j in range(len(lst)):
        for i in range(len(lst[j])):
            if i % 12 == 0 or (i - 4) % 12 == 0:
                pass
            else:
                if lst[j][i] != 'None' and lst[j][i] != '--':
                    lst[j][i] = float(lst[j][i])

    return lst


def write_csv(lst, name='temp.csv'):
    """
    Write csv file from nested list

    :param name:
    :param lst:
    :return:
    """
    with open(name, 'w', newline='') as csv_file:
        writer = csv.writer(csv_file)

        for row in lst:
            row_with_date = None
            writer.writerow(row)


if __name__ == '__main__':
    loc = ''

    filename = 'values_' + loc + '.txt'
    csv_name = 'weather_' + loc + '.csv'

    print('-- Initiating process --')
    l = read_txt(filename)

    print('File read successfully')

    l_int = str_to_num(l)

    print("Values converted to numbers")

    write_csv(l_int, csv_name)

    #farenheit_to_celsius(csv_name)

    print("CSV file created")
