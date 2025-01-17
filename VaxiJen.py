#!/usr/bin/env python3

#Import Libraries
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.select import Select
from selenium.webdriver.chrome.options import Options
import time
import pyperclip
import os
import argparse

# Initialize the parser
parser = argparse.ArgumentParser(description='VaxiJen web scrapping')

# Add arguments
parser.add_argument('chromedriver_path', type=str, help='The path to the chromedriver')
parser.add_argument('output_dir', type=str, help='The output directory')
parser.add_argument('output_file', type=str, help='The output file name')

# Parse the arguments
args = parser.parse_args()

#Navigation options
options =  webdriver.ChromeOptions()
options.add_argument("--no-sandbox")
options.add_argument('--start-maximized')
options.add_argument('--disable-extensions')
options.add_argument('--disable-dev-shm-usage')
options.page_load_strategy = 'none'
options.add_argument("enable-automation")
options.add_argument("--disable-browser-side-navigation")
options.add_argument("--disable-gpu")
options.add_argument("--disable-infobars")     

ser = Service(chromedriver_path)
driver = webdriver.Chrome(service=ser, options=options)

print('Get the target URL')
driver.get('https://www.ddg-pharmfac.net/vaxijen3/')

print('Wait for the webpage to load completely')
time.sleep(30)

print('Find')
Button = driver.find_element(By.NAME, 'uploaded_file')
Button.send_keys(output_dir + output_file)

print('Submit')
driver.find_element(By.NAME, 'submit').click()
time.sleep(200)

print('Copy')
webdriver.ActionChains(driver).key_down(Keys.CONTROL).send_keys("a").perform()
webdriver.ActionChains(driver).key_down(Keys.CONTROL).send_keys("c").perform()

print('Paste')
content = pyperclip.paste()
outFile = open(output_dir + output_file, 'w')
outFile.write(content)
outFile.close()
time.sleep(30)

print('Log out')
driver.quit()

