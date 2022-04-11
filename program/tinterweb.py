#!/usr/bin/python3

import multiprocessing
import os
import urllib.request
import requests


_DOWNLOADS_PATH = "downloads"
_LOCK = multiprocessing.Lock()


def request_json(url, data=None):
    data = data or {}
    with _LOCK:
        print("Requesting:", url)
        if data:
            print("With data:", data)
            response = requests.post(url, data=data)
        else:
            response = requests.get(url)
        return response.json()


def download_file(filename, url):
    path = os.path.join(_DOWNLOADS_PATH, filename)
    with _LOCK:
        if not os.path.exists(path):
            os.makedirs(_DOWNLOADS_PATH, exist_ok=True)
            print("Downloading:", url)
            urllib.request.urlretrieve(url, path)
    return path
