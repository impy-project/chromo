from chromo.util import get_all_models, _cached_data_dir

urls = set()
for Model in get_all_models():
    if hasattr(Model, "_data_url"):
        url = Model._data_url
        urls.add(url)

for url in urls:
    print(f"Downloading/installing from {url}")
    _cached_data_dir(url)
