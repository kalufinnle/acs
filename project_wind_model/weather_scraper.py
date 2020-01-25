# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 14:33:37 2020

@author: Marcus Grunnesjo

Collects weather data from wunderground.com and stores it in a txt file with each row as a day.
Note: to run you must change loc variable in 3 locations to the weather station identifier.


:param loc must be changed

"""

import scrapy
from scrapy.crawler import CrawlerProcess
from scrapy.http import Request
import time


class MySpider(scrapy.Spider):

    name = "wunder_spider"
    allowed_domains = ['wunderground.com']
    loc = 'KCACANTI2'

    start_urls = ['https://www.wunderground.com/dashboard/pws/' + loc + '/table/2019-01-13/2019-01-13/daily']

    def parse(self, response):
        """
        Function calling the parse_page function on the correct pages, by year, month and day
        :param response:
        :return:
        """

        # Base URL
        loc = 'KCACANTI2'
        locator = 'https://www.wunderground.com/dashboard/pws/' + loc + '/table'

        for year in range(2017, 2020):

            for month in range(1, 13):

                days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

                for day in range(days[month - 1] + 1):
                    # Url by date
                    url = locator + '/' + str(year) + '-' + str(month) + '-' + str(day) + '/' + str(year) + '-' + str(
                        month) + '-' + str(day) + '/daily'

                    yield Request(url, callback=self.parse_page)
                    print("Page " + str(year) + '-' + str(month) + '-' + str(day))
                    time.sleep(0.1)

    def parse_page(self, response):
        """
        Function that collects the data values from table on webpage using xpath locators
        :param response:
        :return:
        """

        # Base XPath locator
        loc = 'KCACANTI2'
        locator = './/body//table[@class="history-table desktop-table"]/tbody'
        row = []

        # Iterates through all rows of table
        rows_per_table = 97
        for i in range(1, rows_per_table):
            base = locator + '/tr[' + str(i) + ']'

            time = base + '/td[1]/strong/text()'
            temp = base + '/td[2]/lib-display-unit/span/span/text()'
            dew = base + '/td[3]/lib-display-unit/span/span/text()'
            humid = base + '/td[4]/lib-display-unit/span/span/text()'
            wind = base + '/td[5]/strong/text()'
            speed = base + '/td[6]/lib-display-unit/span/span/text()'
            gust = base + '/td[7]/lib-display-unit/span/span/text()'
            pressure = base + '/td[8]/lib-display-unit/span/span/text()'
            precip_rate = base + '/td[9]/lib-display-unit/span/span/text()'
            precip_acc = base + '/td[10]/lib-display-unit/span/span/text()'
            uv = base + '/td[11]/strong/text()'
            solar = base + '/td[12]/strong/text()'

            # Dictionary of values
            d = {
                'time': response.xpath(time).get(),
                'temp': response.xpath(temp).get(),
                'dew': response.xpath(dew).get(),
                'humid': response.xpath(humid).get(),
                'wind': response.xpath(wind).get(),
                'speed': response.xpath(speed).get(),
                'gust': response.xpath(gust).get(),
                'pressure': response.xpath(pressure).get(),
                'precip_rate': response.xpath(precip_rate).get(),
                'precip_acc': response.xpath(precip_acc).get(),
                'uv': response.xpath(uv).get(),
                'solar': response.xpath(solar).get(),
            }

            # List from dictionary
            row.extend(list(d.values()))

        # Append list to txt file
        filename = 'values_' + loc + '.txt'
        with open(filename, 'a') as file:
            file.write(str(row))
            file.write('\n')


process = CrawlerProcess({
    'USER_AGENT': 'Mozilla/4.0 (compatible; MSIE 7.0; Windows NT 5.1)'
})

process.crawl(MySpider)
process.start()
