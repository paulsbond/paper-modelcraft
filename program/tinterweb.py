#!/usr/bin/python3

import json
import os
import urllib.request
import requests
import multiprocessing


_REQUEST_JSON_CACHE = {}


def request_json(url, data=None):
    if not _REQUEST_JSON_CACHE and os.path.exists("request_cache.json"):
        with open("request_cache.json") as stream:
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
        with open("request_cache.json", "w") as stream:
            json.dump(_REQUEST_JSON_CACHE, stream)
        return response.json()


def download_file(filename, url):
    path = os.path.join("downloads", filename)
    with multiprocessing.Lock():
        if not os.path.exists(path):
            os.makedirs("downloads", exist_ok=True)
            urllib.request.urlretrieve(url, path)
    return path
