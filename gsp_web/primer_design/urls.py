from django.urls import path
from . import views
from .chat_views import chat

urlpatterns = [
    path('chat/', chat, name='chat'),
    path('download-selected/', views.download_selected, name='download_selected'),
    path('kasp/', views.kasp_design, name='kasp_design'),
    path('kasp/download/<str:job_id>/', views.kasp_download_all, name='kasp_download_all'),
    path('kasp/download-selected/', views.kasp_download_selected, name='kasp_download_selected'),
    path('detect-snps/', views.detect_snps, name='detect_snps'),
    path('',              views.index,       name='index'),
    path('run/',          views.run_pipeline, name='run_pipeline'),
    path('progress/<str:job_id>/', views.job_progress, name='job_progress'),
    path('primer-locations/', views.primer_locations, name='primer_locations'),
    path('download/<str:job_id>/', views.download_csv, name='download_csv'),
]
