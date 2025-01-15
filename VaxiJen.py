#!/usr/bin/env python3

print('Import Library')
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.select import Select
from selenium.webdriver.chrome.options import Options
import time
import pyperclip
import os

print('Navigation options')
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

# Modify here!
ser = Service('/path/to/chromedriver')
driver = webdriver.Chrome(service=ser, options=options)

print('Get the target URL')
driver.get('https://www.ddg-pharmfac.net/vaxijen3/')

print('Wait for the webpage to load completely')
time.sleep(30)

print('Find')
Button = driver.find_element(By.NAME, 'uploaded_file')
# Modify also here!
Button.send_keys('/path/to/TempDir/protein.part')

print('Submit')
driver.find_element(By.NAME, 'submit').click()
time.sleep(200)

print('Copy')
webdriver.ActionChains(driver).key_down(Keys.CONTROL).send_keys("a").perform()
webdriver.ActionChains(driver).key_down(Keys.CONTROL).send_keys("c").perform()

print('Paste')
content = pyperclip.paste()
outFile = open('VaxiJen3_predictions.csv', 'w')
outFile.write(content)
outFile.close()
time.sleep(30)

print('Log out')
driver.quit()


