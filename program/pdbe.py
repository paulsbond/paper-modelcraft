import tinterweb


def search(query, filter_list):
    # Documentation: https://www.ebi.ac.uk/pdbe/api/doc/search.html
    url = "https://www.ebi.ac.uk/pdbe/search/pdb/select?"
    data = {"q": query, "fl": filter_list, "rows": 1000000, "wt": "json"}
    response = tinterweb.request_json(url, data)
    return response.get("response", {}).get("docs", [])
