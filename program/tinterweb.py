#!/usr/bin/python3

import json
import multiprocessing
import os
import urllib.request
import requests


_DOWNLOADS_PATH = "downloads"
_REQUEST_JSON_CACHE_PATH = "request_cache.json"
_REQUEST_JSON_CACHE = {}
_LOCK = multiprocessing.Lock()


def request_json(url, data=None):
    data = data or {}
    key = str((url, tuple(sorted(data.items()))))
    with _LOCK:
        if not _REQUEST_JSON_CACHE and os.path.exists(_REQUEST_JSON_CACHE_PATH):
            with open(_REQUEST_JSON_CACHE_PATH) as stream:
                _REQUEST_JSON_CACHE.update(json.load(stream))
        if key in _REQUEST_JSON_CACHE:
            return _REQUEST_JSON_CACHE[key]
        print("Requesting:", url)
        if data:
            print("With data:", data)
            response = requests.post(url, data=data)
        else:
            response = requests.get(url)
        _REQUEST_JSON_CACHE[key] = response.json()
        with open(_REQUEST_JSON_CACHE_PATH, "w") as stream:
            json.dump(_REQUEST_JSON_CACHE, stream)
        return response.json()


def download_file(filename, url):
    path = os.path.join(_DOWNLOADS_PATH, filename)
    with _LOCK:
        if not os.path.exists(path):
            os.makedirs(_DOWNLOADS_PATH, exist_ok=True)
            print("Downloading:", url)
            urllib.request.urlretrieve(url, path)
    return path
