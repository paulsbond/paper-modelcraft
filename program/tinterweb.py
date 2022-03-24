#!/usr/bin/python3

import json
import multiprocessing
import os
import urllib.request
import requests


_DOWNLOADS_PATH = "downloads"
_CACHE_PATH = "request_cache.json"
_LOCK = multiprocessing.Lock()


def request_json(url, data=None):
    data = data or {}
    key = str((url, tuple(sorted(data.items()))))
    with _LOCK:
        cache = {}
        if os.path.exists(_CACHE_PATH):
            with open(_CACHE_PATH) as stream:
                cache = json.load(stream)
        if key in cache:
            return cache[key]
        print("Requesting:", url)
        if data:
            print("With data:", data)
            response = requests.post(url, data=data)
        else:
            response = requests.get(url)
        cache[key] = response.json()
        with open(_CACHE_PATH, "w") as stream:
            json.dump(cache, stream)
        return response.json()


def download_file(filename, url):
    path = os.path.join(_DOWNLOADS_PATH, filename)
    with _LOCK:
        if not os.path.exists(path):
            os.makedirs(_DOWNLOADS_PATH, exist_ok=True)
            print("Downloading:", url)
            urllib.request.urlretrieve(url, path)
    return path
