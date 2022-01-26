#!/usr/bin/python3

import json
import os
import urllib.request
import requests
import multiprocessing


_DOWNLOADS_PATH = "downloads"
_REQUEST_JSON_CACHE_PATH = "request_cache.json"
_REQUEST_JSON_CACHE = {}


def request_json(url, data=None):
    if not _REQUEST_JSON_CACHE and os.path.exists(_REQUEST_JSON_CACHE_PATH):
        with open(_REQUEST_JSON_CACHE_PATH) as stream:
            _REQUEST_JSON_CACHE = json.load(stream)
    data = data or {}
    hash_ = hash((url, tuple(sorted(data.items()))))
    with multiprocessing.Lock():
        if hash_ in _REQUEST_JSON_CACHE:
            return _REQUEST_JSON_CACHE[hash_]
        if data is None:
            response = requests.get(url)
        else:
            response = requests.post(url, data=data)
        _REQUEST_JSON_CACHE[hash_] = response.json()
        with open(_REQUEST_JSON_CACHE_PATH, "w") as stream:
            json.dump(_REQUEST_JSON_CACHE, stream)
        return response.json()


def download_file(filename, url):
    path = os.path.join(_DOWNLOADS_PATH, filename)
    with multiprocessing.Lock():
        if not os.path.exists(path):
            os.makedirs(_DOWNLOADS_PATH, exist_ok=True)
            urllib.request.urlretrieve(url, path)
    return path
